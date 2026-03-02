################################################################################
# Copyright (c) 2026: Quentin Lété
################################################################################

"""
    epec_market_power_model

Economic equilibrium on electricity market for Germany with nodal DC power flow and strategic bidding.
- Upper level: largest company in DE maximizes profit via strategic bids
- Lower level: Market operator clears nodal market via cost-minimizing dispatch

"""
function epec_market_power_model(H)
    dir = joinpath(@__DIR__, "data", "CWE2018_daily")
    buses = CSV.read(joinpath(dir, "Buses.csv"), DataFrame)
    lines = CSV.read(joinpath(dir, "Lines.csv"), DataFrame)
    gens = CSV.read(joinpath(dir, "Generators.csv"), DataFrame)
    gensRE = CSV.read(joinpath(dir, "GeneratorsRE.csv"), DataFrame)
    loads = CSV.read(joinpath(dir, "Loads.csv"), DataFrame)
    load_profiles = CSV.read(joinpath(dir, "LoadProfiles2018.csv"), DataFrame)
    ren_profiles = CSV.read(joinpath(dir, "RenewableProfiles.csv"), DataFrame)
    # Germany filter (Country column assumed)
    buses_DE = filter(r -> r.Country == "D", buses)

    bus_set = Set(string.(buses_DE.Bus))

    lines_DE =
        filter(r -> string(r.FromBus) in bus_set && string(r.ToBus) in bus_set, lines)
    gens_DE = filter(r -> string(r.BusGenerator) in bus_set, gens)
    gensRE_DE = filter(r -> string(r.BusGeneratorRE) in bus_set, gensRE)
    loads_DE = filter(r -> string(r.BusLoad) in bus_set, loads)
    ## ASSIGN COMPANIES TO GENERATORS
    pp = CSV.read(joinpath(dir, "conventional_power_plants_DE.csv"), DataFrame)
    pp.id = string.(pp.id)

    # map generator -> company
    company_map = Dict(pp.id[i] => pp.company[i] for i = 1:nrow(pp))

    default_company = "OTHER"
    gens_DE.company =
        [get(company_map, string(g), default_company) for g in gens_DE.Generator]

    ## IDENTIFY LARGEST COMPANY BY CAPACITY
    cap_by_company = combine(groupby(gens_DE, :company), :Capacity => sum => :TotalCapacity)
    sort!(cap_by_company, :TotalCapacity, rev = true)

    largest_company = cap_by_company.company[1]
    println("Largest company in Germany by capacity: ", largest_company)

    # index sets
    bus_ids = collect(string.(buses_DE.Bus))
    bus_to_idx = Dict(b => i for (i, b) in enumerate(bus_ids))

    ng = nrow(gens_DE)
    nb = length(bus_ids)
    nl = nrow(lines_DE)

    # ---------------- Battery data extraction ----------------
    # identify batteries by technology (explicit: "Storage technologies")
    batt_pp = filter(r -> occursin("storage technologies", lowercase(string(r.technology))), pp)

    # use net capacity as power capacity (MW)
    batt_P = Float64.(batt_pp.capacity_net_bnetza)
    batt_E = batt_P .* 1.0 # E/P = 1
    eta_rt = 0.9
    eta_c = sqrt(eta_rt)
    eta_d = sqrt(eta_rt)

    # ---------------- Assign batteries to closest bus ----------------
    bus_lat = Dict(string(buses_DE.Bus[i]) => buses_DE.lat_pred[i] for i = 1:nrow(buses_DE))
    bus_lon = Dict(string(buses_DE.Bus[i]) => buses_DE.lon_pred[i] for i = 1:nrow(buses_DE))
    function haversine(lat1, lon1, lat2, lon2)
        R = 6371.0
        dlat = deg2rad(lat2 - lat1)
        dlon = deg2rad(lon2 - lon1)
        a = sin(dlat/2)^2 + cos(deg2rad(lat1))*cos(deg2rad(lat2))*sin(dlon/2)^2
        return 2R * asin(sqrt(a))
    end

    batt_bus = Int[]
    for r in eachrow(batt_pp)
        dmin = Inf;
        best = 0
        for (b, latb) in bus_lat
            d = haversine(r.lat, r.lon, latb, bus_lon[b])
            if d < dmin
                dmin = d;
                best = bus_to_idx[b]
            end
        end
        push!(batt_bus, best)
    end

    nbatt = length(batt_bus)

    hours = 1:H

    ## DEFINE BILEVEL MODEL
    model = BilevelModel(Ipopt.Optimizer, mode = BilevelJuMP.ComplementMode())

    # ---------------- Upper level variables (strategic firm bids) ----------------
    strategic_gens = findall(gens_DE.company .== largest_company)

    @variable(Upper(model), bid_price[i in strategic_gens, h in hours] >= 0)
    @variable(Upper(model), bid_quantity[i in strategic_gens, h in hours] >= 0)

    # ---------------- Lower level variables ----------------
    @variable(Lower(model), pg[1:ng, hours] >= 0)
    @variable(Lower(model), pr[1:nrow(gensRE_DE), hours] >= 0)
    @variable(Lower(model), pch[1:nbatt, hours] >= 0)
    @variable(Lower(model), pdis[1:nbatt, hours] >= 0)
    @variable(Lower(model), soc[1:nbatt, hours] >= 0)
    @variable(Lower(model), theta[1:nb, hours])
    @variable(Lower(model), f[1:nl, hours])

    # ---------------- Lower level constraints ----------------
    ref_bus = 1
    for h in hours
        @constraint(Lower(model), theta[ref_bus, h] == 0)
    end

    # generator limits
    for i = 1:ng, h in hours
        cap = gens_DE.Capacity[i]
        if i in strategic_gens
            @constraint(Lower(model), pg[i, h] <= bid_quantity[i, h])
            @constraint(Lower(model), bid_quantity[i, h] <= cap) # you can withhold capacity but not exceed it
        else
            @constraint(Lower(model), pg[i, h] <= cap)
        end
    end

    # renewable availability (hour 1, illustrative)
    for (k, r) in enumerate(eachrow(gensRE_DE)), h in hours
        prof = r.DynamicProfileRE
        frac = coalesce(r.FractionREProfile, 1.0)
        avail = frac * ren_profiles[h, prof]
        @constraint(Lower(model), pr[k, h] <= max(avail, 0.0))
    end

    # line data
    from_bus = [bus_to_idx[string(lines_DE.FromBus[i])] for i = 1:nl]
    to_bus = [bus_to_idx[string(lines_DE.ToBus[i])] for i = 1:nl]
    x = Float64.(lines_DE.Reactance)
    limit = Float64.(coalesce.(lines_DE.FlowLimitForw, 1e9))

    for l = 1:nl, h in hours
        @constraint(
            Lower(model),
            f[l, h] == (theta[from_bus[l], h] - theta[to_bus[l], h]) / x[l]
        )
        @constraint(Lower(model), f[l, h] <= limit[l])
        @constraint(Lower(model), -limit[l] <= f[l, h])
    end

    # ---------------- Battery constraints ----------------
    for b = 1:nbatt, h in hours
        @constraint(Lower(model), pch[b, h] <= batt_P[b])
        @constraint(Lower(model), pdis[b, h] <= batt_P[b])
        @constraint(Lower(model), soc[b, h] <= batt_E[b])
    end

    for b = 1:nbatt
        @constraint(Lower(model), soc[b, 1] == eta_c * pch[b, 1] - pdis[b, 1] / eta_d)
        for h = 2:H
            @constraint(
                Lower(model),
                soc[b, h] == soc[b, h-1] + eta_c * pch[b, h] - pdis[b, h] / eta_d
            )
        end
    end

    # build load per bus
    node_load = zeros(nb, length(hours))
    for r in eachrow(loads_DE), h in hours
        b = bus_to_idx[string(r.BusLoad)]
        prof = r.LoadProfile
        frac = r.FractionLoadProfile
        node_load[b, h] += frac * load_profiles[h, prof]
    end

    for b = 1:nb, h in hours
        gens = collect(i for i = 1:ng if bus_to_idx[string(gens_DE.BusGenerator[i])] == b)
        conv = isempty(gens) ? 0.0 : sum(pg[i, h] for i in gens)
        gensRE = collect(
            k for k = 1:nrow(gensRE_DE) if
            bus_to_idx[string(gensRE_DE.BusGeneratorRE[k])] == b
        )
        ren = isempty(gensRE) ? 0.0 : sum(pr[k, h] for k in gensRE)
        batteries = collect(j for j = 1:nbatt if batt_bus[j] == b)
        bat = isempty(batteries) ? 0.0 : sum(pdis[j, h] - pch[j, h] for j in batteries)
        flow = sum(
            (from_bus[l] == b ? f[l, h] : 0.0) - (to_bus[l] == b ? f[l, h] : 0.0) for
            l = 1:nl
        )
        @constraint(Lower(model), conv + ren + bat - node_load[b, h] - flow == 0)
    end

    # ---------------- Lower level objective ----------------
    @objective(
        Lower(model),
        Min,
        sum(
            (i in strategic_gens ? bid_price[i, h] : gens_DE.VariableCost[i]) * pg[i, h] for
            i = 1:ng, h in hours
        )
    )

    # ---------------- Upper level objective (profit) ----------------
    # Profit = revenue - cost
    @objective(
        Upper(model),
        Max,
        sum(
            (bid_price[i, h] - gens_DE.VariableCost[i]) * pg[i, h] for
            i in strategic_gens, h in hours
        )
    )

    return model
end
