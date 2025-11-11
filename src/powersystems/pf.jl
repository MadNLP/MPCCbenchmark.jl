
#=
    Implement powerflow with PV/PQ switches.
=#

"""
    poweflow_model(data::PowerData)

Implement powerflow with PV/PQ switches from [1].

# Reference

[1] Murray, W., Tinoco De Rubira, T.,  Wigington, A.
"A robust and informative method for solving large-scale power flow problems."
Computational optimization and applications, 62(2), 431-475, 2015.

"""
function powerflow_model(data; load_factor=1.0)
    nbus = length(data.bus)
    ngen = length(data.gen)
    nbranch = length(data.branch)
    narc = length(data.arc)

    pq_buses = [k for k in 1:nbus if data.bus[k].type == 1]
    pv_buses = [k for k in 1:nbus if data.bus[k].type == 2]
    ref_buses = [k for k in 1:nbus if data.bus[k].type == 3]

    # Active power generation
    pg = [g.pg for g in data.gen]

    connected_arcs = [Int[] for k in 1:nbus]
    for (k, arc) in enumerate(data.arc)
        b = arc.bus
        push!(connected_arcs[b], k)
    end
    connected_gens = [Int[] for k in 1:nbus]
    qmin = zeros(nbus)
    qmax = zeros(nbus)
    qg0 = zeros(nbus)
    vm0 = zeros(nbus)
    for (k, gen) in enumerate(data.gen)
        b = gen.bus
        push!(connected_gens[b], k)
        # Aggregate reactive power generations
        qmin[b] += gen.qmin
        qmax[b] += gen.qmax
        qg0[b] += gen.qg
    end
    for b in 1:nbus
        vm0[b] = data.bus[b].vm
    end

    nref, npv, npq = length(ref_buses), length(pv_buses), length(pq_buses)

    model = JuMP.Model()
    JuMP.set_name(model, "Powerflow-with-switches")

    JuMP.@variable(model, va[i in 1:nbus], start=deg2rad(data.bus[i].va))
    JuMP.@variable(model, vm[i in 1:nbus], start=vm0[i])
    JuMP.@variable(model, qmin[i] <= qg[i in pv_buses] <= qmax[i], start=qg0[i])
    JuMP.@variable(model, p[i in 1:narc])
    JuMP.@variable(model, q[i in 1:narc])

    # Complementarity variables for PV/PQ switches
    JuMP.@variable(model, Δ[i in pv_buses])

    # Fix voltage angle and magnitude at ref node
    for b in ref_buses
        JuMP.@constraint(model, va[b] == data.bus[b].va)
        JuMP.@constraint(model, vm[b] == vm0[b])
    end
    # Encode PV/PQ switches at at PV nodes with complementarity constraints
    for b in pv_buses
        @constraint(model, vm[b] == vm0[b]+ Δ[b])
        @constraint(model, [Δ[b], qg[b]] ∈ MOI.Complements(2))
    end

    # Power balance at PQ buses
    for b in pq_buses
        JuMP.@constraint(model,
            sum(p[a] for a in connected_arcs[b]) ==
            - load_factor * data.bus[b].pd
            - data.bus[b].gs * vm[b]^2
        )
        JuMP.@constraint(model,
            sum(q[a] for a in connected_arcs[b]) ==
            - load_factor * data.bus[b].qd
            + data.bus[b].bs * vm[b]^2
        )
    end
    # Power balance at PV buses
    for b in pv_buses
        JuMP.@constraint(model,
            sum(p[a] for a in connected_arcs[b]) ==
            sum(pg[g] for g in connected_gens[b])
            - load_factor * data.bus[b].pd
            - data.bus[b].gs * vm[b]^2
        )
        JuMP.@constraint(model,
            sum(q[a] for a in connected_arcs[b]) ==
            qg[b]
            - load_factor * data.bus[b].qd
            + data.bus[b].bs * vm[b]^2
        )
    end

    # Power flow constraints
    for l in data.branch
        f_idx = l.f_idx
        t_idx = l.t_idx

        p_fr = p[f_idx]
        q_fr = q[f_idx]
        p_to = p[t_idx]
        q_to = q[t_idx]

        vm_fr = vm[l.f_bus]
        vm_to = vm[l.t_bus]
        va_fr = va[l.f_bus]
        va_to = va[l.t_bus]

        # FR
        JuMP.@constraint(
            model,
            p_fr - l.c5*vm_fr^2
            - l.c3 * (vm_fr * vm_to * cos(va_fr - va_to)) - l.c4 * (vm_fr * vm_to * sin(va_fr - va_to))
            == 0.0
        )
        JuMP.@constraint(
            model,
            q_fr + l.c6*vm_fr^2
            + l.c4 * (vm_fr * vm_to * cos(va_fr - va_to))
            - l.c3 * (vm_fr * vm_to * sin(va_fr - va_to))
            == 0.0
        )

        # TO
        JuMP.@constraint(
            model,
            p_to - l.c7*vm_to^2
            - l.c1 * (vm_to * vm_fr * cos(va_to - va_fr))
            - l.c2 * (vm_to * vm_fr * sin(va_to - va_fr))
            == 0.0
        )
        JuMP.@constraint(
            model,
            q_to + l.c8*vm_to^2
            + l.c2 * (vm_to * vm_fr * cos(va_to - va_fr))
            - l.c1 * (vm_to * vm_fr * sin(va_to - va_fr))
            == 0.0
        )
    end

    return model
end

