
#=
    Implement corrective SCOPF.
=#

function _calc_branch_t(branch::Dict{String,<:Any})
    tap_ratio = branch["tap"]
    angle_shift = branch["shift"]

    tr = tap_ratio .* cos.(angle_shift)
    ti = tap_ratio .* sin.(angle_shift)

    return tr, ti
end

function _calc_branch_y(branch::Dict{String,<:Any})
    y = LinearAlgebra.pinv(branch["br_r"] + im * branch["br_x"])
    g, b = real(y), imag(y)
    return g, b
end

"""
    scopf_model(data, contingencies)

Implement corrective AC-SCOPF with complementarity constraints from [1].

# Reference

[1] F Pacaud, A Nurkanović, A Pozharskiy, A Montoison, S Shin.
"An Augmented Lagrangian Method on GPU for Security-Constrained AC Optimal Power Flow", 2025.

"""
function scopf_model(
    data,
    contingencies;
    scale_obj=1e-4,
    load_factor=1.0,
)
    ngen = length(data[:gen])
    alpha = ones(ngen)

    # Parse contingencies
    K = length(contingencies) + 1

    # Build model
    model = JuMP.Model()
    JuMP.set_name(model, "Corrective-AC-SCOPF")

    JuMP.@variable(model, va[i in keys(data[:bus]), 1:K])
    JuMP.@variable(model, data[:bus][i]["vmin"] <= vm[i in keys(data[:bus]), 1:K] <= data[:bus][i]["vmax"], start=1.0)
    JuMP.@variable(model, data[:gen][i]["pmin"] <= pg[i in keys(data[:gen]), 1:K] <= data[:gen][i]["pmax"])
    JuMP.@variable(model, data[:gen][i]["qmin"] <= qg[i in keys(data[:gen]), 1:K] <= data[:gen][i]["qmax"])
    JuMP.@variable(model, -data[:branch][l]["rate_a"] <= p[(l,i,j) in data[:arcs], 1:K] <= data[:branch][l]["rate_a"])
    JuMP.@variable(model, -data[:branch][l]["rate_a"] <= q[(l,i,j) in data[:arcs], 1:K] <= data[:branch][l]["rate_a"])
    JuMP.@variable(model, 0.0 <= ρp[i in keys(data[:gen]), 2:K])
    JuMP.@variable(model, 0.0 <= ρn[i in keys(data[:gen]), 2:K])
    JuMP.@variable(model, 0.0 <= vp[i in keys(data[:gen]), 2:K])
    JuMP.@variable(model, 0.0 <= vn[i in keys(data[:gen]), 2:K])

    # Automatic adjustment of generators
    JuMP.@variable(model, Δ[1:K-1])

    JuMP.@objective(model, Min, scale_obj * sum(gen["cost"][1]*pg[i, 1]^2 + gen["cost"][2]*pg[i, 1] + gen["cost"][3] for (i,gen) in data[:gen]))

    for k in 2:K
        # Droop control
        for (j, i) in enumerate(keys(data[:gen]))
            pmin, pmax = data[:gen][i]["pmin"], data[:gen][i]["pmax"]
            @constraint(model, ρp[i, k] - ρn[i, k] == pg[i, k] - pg[i, 1] - alpha[j] * Δ[k-1])
            @constraint(model, [(pmax - pg[i, k]), ρn[i, k]] in MOI.Complements(2))
            @constraint(model, [(pg[i, k] - pmin), ρp[i, k]] in MOI.Complements(2))
        end
        # Voltage magnitude are not adjusted at PV buses
        for g in keys(data[:gen])
            b = data[:gen][g]["gen_bus"]
            qmin, qmax = data[:gen][g]["qmin"], data[:gen][g]["qmax"]
            @constraint(model, vp[g, k] - vn[g, k] == vm[b, k] - vm[b, 1])
            if isfinite(qmax)
                @constraint(model, [(qmax - qg[g, k]), vn[g, k]] in MOI.Complements(2))
            end
            if isfinite(qmin)
                @constraint(model, [(qg[g, k] - qmin), vp[g, k]] in MOI.Complements(2))
            end
        end
        # Set flux to 0
        for (l, i, j) in data[:arcs]
            if l == contingencies[k-1]
                @constraint(model, p[(l, i, j), k] == 0.0)
                @constraint(model, q[(l, i, j), k] == 0.0)
            end
        end
    end

    for k in 1:K
        for (i, bus) in data[:ref_buses]
            JuMP.@constraint(model, va[i, k] == 0)
        end

        for (i,bus) in data[:bus]
            bus_loads = [data[:load][l] for l in data[:bus_loads][i]]
            bus_shunts = [data[:shunt][s] for s in data[:bus_shunts][i]]

            JuMP.@constraint(model,
                sum(p[a, k] for a in data[:bus_arcs][i]) ==
                sum(pg[g, k] for g in data[:bus_gens][i]) -
                sum(load_factor * load["pd"] for load in bus_loads) -
                sum(shunt["gs"] for shunt in bus_shunts)*vm[i, k]^2
            )

            JuMP.@constraint(model,
                sum(q[a, k] for a in data[:bus_arcs][i]) ==
                sum(qg[g, k] for g in data[:bus_gens][i]) -
                sum(load_factor * load["qd"] for load in bus_loads) +
                sum(shunt["bs"] for shunt in bus_shunts)*vm[i, k]^2
            )
        end

        # Branch power flow physics and limit constraints
        for (i,branch) in data[:branch]
            if (k >= 2) && i == contingencies[k-1]
                continue
            end
            f_idx = (i, branch["f_bus"], branch["t_bus"])
            t_idx = (i, branch["t_bus"], branch["f_bus"])

            p_fr = p[f_idx, k]
            q_fr = q[f_idx, k]
            p_to = p[t_idx, k]
            q_to = q[t_idx, k]

            vm_fr = vm[branch["f_bus"], k]
            vm_to = vm[branch["t_bus"], k]
            va_fr = va[branch["f_bus"], k]
            va_to = va[branch["t_bus"], k]

            g, b = _calc_branch_y(branch)
            tr, ti = _calc_branch_t(branch)
            ttm = tr^2 + ti^2
            g_fr = branch["g_fr"]
            b_fr = branch["b_fr"]
            g_to = branch["g_to"]
            b_to = branch["b_to"]

            # From side of the branch flow
            JuMP.@constraint(model, p_fr ==  (g+g_fr)/ttm*vm_fr^2 + (-g*tr+b*ti)/ttm*(vm_fr*vm_to*cos(va_fr-va_to)) + (-b*tr-g*ti)/ttm*(vm_fr*vm_to*sin(va_fr-va_to)) )
            JuMP.@constraint(model, q_fr == -(b+b_fr)/ttm*vm_fr^2 - (-b*tr-g*ti)/ttm*(vm_fr*vm_to*cos(va_fr-va_to)) + (-g*tr+b*ti)/ttm*(vm_fr*vm_to*sin(va_fr-va_to)) )

            # To side of the branch flow
            JuMP.@constraint(model, p_to ==  (g+g_to)*vm_to^2 + (-g*tr-b*ti)/ttm*(vm_to*vm_fr*cos(va_to-va_fr)) + (-b*tr+g*ti)/ttm*(vm_to*vm_fr*sin(va_to-va_fr)) )
            JuMP.@constraint(model, q_to == -(b+b_to)*vm_to^2 - (-b*tr+g*ti)/ttm*(vm_to*vm_fr*cos(va_to-va_fr)) + (-g*tr-b*ti)/ttm*(vm_to*vm_fr*sin(va_to-va_fr)) )

            # Apparent power limit, from side and to side
            JuMP.@constraint(model, p_fr^2 + q_fr^2 <= branch["rate_a"]^2)
            JuMP.@constraint(model, p_to^2 + q_to^2 <= branch["rate_a"]^2)
        end
    end

    return model
end
