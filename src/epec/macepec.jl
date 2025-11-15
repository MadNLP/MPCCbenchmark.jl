#=
    A set of EPEC instances from MacEPEC

    https://wiki.mcs.anl.gov/leyffer/index.php/MacEPEC

    JuMP translation with ChatGPT.
=#


#######################################################################
# MPEC of example due to Fukushima and Pang, "Quasi-Variational
# Inequalities, Nash-Equilibria, and Multi-Leader-Follower Games".
#
# The generator (firm) is the Stackelberg leader, and the ISO,
# arbitrager, and market clearing are the followers in this game.
# Another possibility is to have the ISO as a Leader as well.
#
# ampl coding by S. Leyffer, Argonne National Laboratory, Jan. 2005.
#######################################################################
struct MacEPECElectricityMarket
    # nodes
    N::UnitRange{Int}
    # firms
    F::UnitRange{Int}
    ARCS::Vector{Tuple{Int, Int}}
    # parameters
    e::Dict{Tuple{Int, Int}, Float64}
    # Generation cost c and capacity CAP by (f,i)
    c::Dict{Tuple{Int, Int}, Float64}
    CAP::Dict{Tuple{Int, Int}, Float64}
    # Price and quantity intercepts
    P0::Dict{Int, Float64}
    Q0::Dict{Int, Float64}
end

function _electricity_model(data::MacEPECElectricityMarket)
    ARCS = data.ARCS
    F = data.F
    N = data.N
    e = data.e
    P0, Q0 = data.P0, data.Q0
    m = Model()

    # amount produced by f at node n1, sold at n2
    @variable(m, s[f in F, (i,j) in ARCS] >= 0)
    # amount of shipment from n1 to n2
    @variable(m, y[(i,j) in ARCS] >= 0)
    # total sales at node n
    @variable(m, S[i in N] >= 0)
    # slacks for easier complementarity
    @variable(m, ss[(i,j) in ARCS] >= 0)
    # amount bought by arbitrager at n1, sold at n2
    @variable(m, a[(i,j) in ARCS] >= 0)
    # unit charge of shipment received by ISO
    @variable(m, w[(i,j) in ARCS], start=1.0)
    @variable(m, sw[(i,j) in ARCS] >= 0.0)

    # multipliers (each firm has different multipliers)
    @variable(m, l_cap[f in F, i in N] >= 0)           # capacity constraint
    @variable(m, l_sal[f in F, i in N])                # total sales
    @variable(m, l_slk[f in F, (i,j) in ARCS])         # slack multipliers (free)
    @variable(m, l_arb[f in F, (i,j) in ARCS] >= 0)    # arbitrager multiplier >=0
    @variable(m, l_mar[f in F, (i,j) in ARCS])         # market clearing multiplier (free)
    @variable(m, l_ISO[f in F, (i,j) in ARCS] >= 0)    # ISO multiplier >=0

    # --- objective: minimize MultNorm = sum l_arb + sum l_ISO (same as AMPL)
    @objective(m, Min, sum(l_arb[f,(i,j)] for f in F, (i,j) in ARCS) + sum(l_ISO[f,(i,j)] for f in F, (i,j) in ARCS))

    ### first order equations wrt s[f,i,j]
    @expression(m, comp1[f in F, (i,j) in ARCS], - (P0[j] - P0[j]/Q0[j] * S[j]) + w[(i,j)] + data.c[f,i] + l_cap[f,i] - l_sal[f,j] - l_mar[f,(i,j)])
    @constraint(m, FO_s[f in F, (i,j) in ARCS], [comp1[f, (i,j)], s[f, (i,j)]] in MOI.Complements(2))

    ### first order equations wrt y[i,j] for every leader f
    @expression(m, comp2[f in F, (i,j) in ARCS], l_mar[f,(i,j)] + l_ISO[f,(i,j)] * (-w[(i,j)] + e[i,j]))
    @constraint(m, FO_y[f in F, (i,j) in ARCS], [comp2[f, (i,j)], y[(i,j)]] in MOI.Complements(2))

    ### first order equations wrt S[i] for every leader f
    @expression(
        m,
        comp3[f in F, i in N],
        P0[i]/Q0[i] * sum(s[f, (j,i)] + a[(j,i)] - a[(i, j)] for (j, ii) in ARCS if ii == i) +
        l_sal[f,i] +
        P0[i]/Q0[i]*sum(l_slk[f, (j,i)] for (j, ii) in ARCS if ii == i) +
        P0[i]/Q0[i]*sum(l_slk[f, (i,j)] for (ii, j) in ARCS if ii == i)
    )
    @constraint(m, FO_S[f in F, i in N], [comp3[f, i], S[i]] in MOI.Complements(2))

    ### first order equations wrt ss[i,j] for every leader f
    @expression(m, comp4[f in F, (i,j) in ARCS], l_slk[f,(i,j)] + l_arb[f,(i,j)] * a[(i,j)])
    @constraint(m, FO_ss[f in F, (i, j) in ARCS], [comp4[f, (i,j)], ss[(i, j)]] in MOI.Complements(2))

    ### first order equations wrt a[i,j] for every leader f
    @expression(m, comp5[f in F, (i,j) in ARCS], -(P0[j] - P0[j]/Q0[j]*S[j]) + (P0[i] - P0[i]/Q0[i]*S[i]) + w[(i, j)] - l_sal[f, j] + l_sal[f, i] + l_arb[f, (i, j)] * ss[(i, j)] - l_mar[f, (i, j)])
    @constraint(m, FO_a[f in F, (i, j) in ARCS], [comp5[f, (i,j)], a[(i,j)]] in MOI.Complements(2))

    ### first order equations wrt w[i,j] for every leader f
    @expression(m, comp6[f in F, (i,j) in ARCS], (a[(i,j)] + s[f,(i,j)]) - l_slk[f,(i,j)] - l_ISO[f,(i,j)]*y[(i,j)])
    @constraint(m, FO_w[f in F, (i, j) in ARCS], [comp6[f, (i, j)], sw[(i, j)]] in MOI.Complements(2))

    ### capacity constraint
    @expression(m, comp7[f in F, i in N], data.CAP[f, i] - sum(s[f, (i, j)] for (ii, j) in ARCS if ii == i))
    @constraint(m, [f in F, i in N], [comp7[f, i], l_cap[f, i]] in MOI.Complements(2))

    ### total sales at node i
    @constraint(m, sales[i in N], -S[i] + sum(s[f, (j,i)] for f in F, (j, ii) in ARCS if ii == i) + sum(a[(j,i)]-a[(i,j)] for (ii, j) in ARCS if ii == i) == 0.0)

    ### define slacks for complementarity
    @constraint(m, slacks[(i, j) in ARCS], - ss[(i,j)] - (P0[j] - P0[j]/Q0[j]*S[j]) + (P0[i] - P0[i]/Q0[i]*S[i]) + w[(i,j)] == 0.0)
    @constraint(m, swcns[(i, j) in ARCS], sw[(i, j)] == e[(i, j)] - w[(i, j)])

    ### market clearing condition
    @constraint(m, [(i, j) in ARCS], y[(i, j)] - sum(s[f, (i, j)] for f in F) - a[(i, j)] == 0.0)

    ### arbitrager's optimality conditions (follower)
    @constraint(m, [(i, j) in ARCS], [a[(i, j)], ss[(i, j)]] in MOI.Complements(2))

    ### ISO's optimality conditions
    @constraint(m, [(i, j) in ARCS], [y[(i, j)], sw[(i, j)]] in MOI.Complements(2))

    return m
end

function macepec_electric002_model()
    data = MacEPECElectricityMarket(
        1:3,
        1:2,
        [(1,1), (1,2), (1,3), (2,1), (2,2), (2,3), (3,1), (3,2), (3,3)],
        Dict(
            (1,1)=>0.0,  (1,2)=>5.0,  (1,3)=>6.0,
            (2,1)=>4.5,  (2,2)=>0.0,  (2,3)=>3.0,
            (3,1)=>4.1,  (3,2)=>2.3,  (3,3)=>0.0
        ),
        Dict(
            (1,1)=>15, (1,2)=>15, (1,3)=>15,
            (2,1)=>15, (2,2)=>15, (2,3)=>15
        ),
        Dict(
            (1,1)=>100, (1,2)=>50,  (1,3)=>0,
            (2,1)=>0,   (2,2)=>100, (2,3)=>50
        ),
        Dict(1=>40, 2=>35, 3=>32),
        Dict(1=>500, 2=>400, 3=>600),
    )
    model = _electricity_model(data)
    JuMP.set_name(model, "MacEpec-electric002")
    return model
end

function macepec_electric004_model()
    data = MacEPECElectricityMarket(
        1:14,
        1:2,
        [
            (1,13), (1,14), (2,3), (2,5), (2,10), (3,2), (3,5), (3,14), (4,13),
            (5,2), (5,3), (5,7), (5,11), (6,9), (6,10), (7,5), (7,10), (7,11),
            (8,11), (8,12), (9,6), (9,12), (10,2), (10,6), (10,7), (10,13), (11,5),
            (11,7), (11,8), (12,8), (12,9), (13,1), (13,4), (13,10), (14,1), (14,3)
        ],
        Dict(
            (1,13)=>8.2274, (1,14)=>6.3038,
            (2,3)=>7.0075, (2,5)=>6.2171, (2,10)=>7.5869,
            (3,2)=>7.0075, (3,5)=>6.7488, (3,14)=>7.4801,
            (4,13)=>7.6682,
            (5,2)=>6.2171, (5,3)=>6.7488, (5,7)=>6.9217, (5,11)=>7.8892,
            (6,9)=>9.3093, (6,10)=>7.5246,
            (7,5)=>6.9217, (7,10)=>6.9882, (7,11)=>7.9131,
            (8,11)=>7.0559, (8,12)=>5.8929,
            (9,6)=>9.3093, (9,12)=>7.4855,
            (10,2)=>7.5869, (10,6)=>7.5246, (10,7)=>6.9882, (10,13)=>6.9950,
            (11,5)=>7.8892, (11,7)=>7.9131, (11,8)=>7.0559,
            (12,8)=>5.8929, (12,9)=>7.4855,
            (13,1)=>8.2274, (13,4)=>7.6682, (13,10)=>6.9950,
            (14,1)=>6.3038, (14,3)=>7.4801
        ),
        Dict(
            (1,1)=>15, (1,2)=>16, (1,3)=>13, (1,4)=>14, (1,5)=>17, (1,6)=>16, (1,7)=>15,
            (1,8)=>0,  (1,9)=>0,  (1,10)=>0, (1,11)=>0, (1,12)=>0, (1,13)=>0, (1,14)=>0,
            (2,1)=>0,  (2,2)=>0,  (2,3)=>0,  (2,4)=>0,  (2,5)=>0,  (2,6)=>0,  (2,7)=>0,
            (2,8)=>14, (2,9)=>15, (2,10)=>12, (2,11)=>11, (2,12)=>15, (2,13)=>13, (2,14)=>14
        ),
        Dict(
            (1,1)=>2500, (1,2)=>4200, (1,3)=>2990, (1,4)=>4000, (1,5)=>3000,
            (1,6)=>4100, (1,7)=>2700, (1,8)=>0, (1,9)=>0, (1,10)=>0, (1,11)=>0, (1,12)=>0, (1,13)=>0, (1,14)=>0,
            (2,1)=>0, (2,2)=>0, (2,3)=>0, (2,4)=>0, (2,5)=>0, (2,6)=>0, (2,7)=>0,
            (2,8)=>2000, (2,9)=>5500, (2,10)=>2900, (2,11)=>2100, (2,12)=>5000, (2,13)=>3100, (2,14)=>5100
        ),
        Dict(
            1=>413.26, 2=>385.93, 3=>364.48, 4=>431.31, 5=>342.84, 6=>429.26, 7=>502.25,
            8=>429.98, 9=>429.32, 10=>429.03, 11=>421.48, 12=>429.92, 13=>345.39, 14=>343.87
        ),
        Dict(
            1=>678.60, 2=>2035.80, 3=>1469.52, 4=>1952.90, 5=>734.76, 6=>1826.80, 7=>1385.60,
            8=>1040.76, 9=>1577.60, 10=>2143.32, 11=>2189.12, 12=>1350.88, 13=>2753.10, 14=>1404.00
        )
    )
    model = _electricity_model(data)
    JuMP.set_name(model, "MacEpec-electric004")
    return model
end

#######################################################################
# ex-001-MPEC.mod: MPEC formulation of EPEC with 2 leaders
#
# leader i's problem is:
#
#     minimize    (x_i+1)^2
#     subject to  s = x_1 + x_2 + y
#                 0 <= s  _|_  y >= 0
#######################################################################
function macepec_ex001_model()
    I = 1:2
    model = Model()
    JuMP.set_name(model, "MacEPEC-ex001")

    @variable(model, x[i in I])
    @variable(model, s >= 0)
    @variable(model, y >= 0)
    @variable(model, lambda[i in I])
    @variable(model, sigma[i in I] >= 0)
    @variable(model, nu[i in I] >= 0)
    @variable(model, xi[i in I] >= 0)

    @objective(model, Min, sum(xi[i] for i in I))
    # First-order conditions
    @constraint(model, [i in I], 2 * (x[i] + 1) + lambda[i] == 0)
    @constraint(model, [i in I], lambda[i] - nu[i] + xi[i] * s == 0)
    @constraint(model, [i in I], -lambda[i] - sigma[i] + xi[i] * y == 0)

    @constraint(model, [i in I],
        [s - (x[1] + x[2] + y), lambda[i]] in MOI.Complements(2))
    @constraint(model, [i in I],
        [s, sigma[i]] in MOI.Complements(2))
    @constraint(model, [i in I],
        [y, nu[i]] in MOI.Complements(2))
    @constraint(model, [s, y] in MOI.Complements(2))

    return model
end

#######################################################################
# ex-4-MPEC.mod
#
# MPEC formulation of Multi-Leader-Follower Game
#
# Example 4 from Fukushima and Pang, "Quasi-Variational Inequalities,
# Generalized Nash Equlibria, and Multi-Leader-Follower Games", to
# appear in Computational Management Science.
#
#######################################################################
function macepec_ex4_model()
    I = 1:2
    model = Model()
    JuMP.set_name(model, "MacEPEC-ex4")
    @variable(model, 0 <= x[i in I] <= 1)
    @variable(model, y >= 0)
    @variable(model, s >= 0)

    @variable(model, chi[i in I])
    @variable(model, psi[i in I] >= 0)
    @variable(model, sigma[i in I] >= 0)
    @variable(model, xi[i in I] >= 0)
    @variable(model, mu[i in I])

    @objective(model, Min, sum(xi[i] for i in I))

    # First-order conditions (leaders)
    @constraint(model, 0.5 - mu[1] == chi[1])
    @constraint(model, 1.0 - mu[1] == psi[1] - xi[1] * s)

    @constraint(model, -0.5 - mu[2] == chi[2])
    @constraint(model, -1.0 - mu[2] == psi[2] - xi[2] * s)

    # First-order conditions wrt slacks s
    @constraint(model, [i in I], mu[i] - sigma[i] + xi[i] * y == 0)

    # Definition of slacks
    @constraint(model, s == -1 + x[1] + x[2] + y)

    # Complementarity conditions per player
    @constraint(model, [i in I],
        [chi[i], x[i]] in MOI.Complements(2))
    @constraint(model, [i in I],
        [y, psi[i]] in MOI.Complements(2))
    @constraint(model, [i in I],
        [s, sigma[i]] in MOI.Complements(2))
    @constraint(model, [s, y] in MOI.Complements(2))

    return model
end

#########################################################################
# outrata3-MPEC.mod
# Original AMPL coding by Sven Leyffer, Argonne National Laboratory, 2004.
#
# MPEC formulation of ...
#
# Multi-Leader-Follower Game (MLF)) derived from outrata31-outrata34,
# see J. Outrata, SIAM J. Optim. 4(2), pp.340ff, 1994. All Stackelberg
# players have the same constraints, but different objectives.
#
#########################################################################
function macepec_outrata3_model()
    nl = 4
    L = 1:nl
    I = 1:4

    model = Model()
    JuMP.set_name(model, "MacEPEC-outrata3")

    # === Leader variables ===
    @variable(model, 0 <= x[l in L] <= 10)
    @expression(model, xa, sum(x[l] for l in L) / nl)
    # === Follower variables ===
    @variable(model, y[i in I] >= 0)
    @variable(model, s[i in I] >= 0)
    # === Multipliers ===
    @variable(model, mu[i in I, l in L])
    @variable(model, xi[i in I, l in L] >= 0)
    @variable(model, sigma[i in I, l in L] >= 0)
    @variable(model, psi[i in I, l in L] >= 0)

    # === Objective ===
    @objective(model, Min, sum(xi[i, l] for i in I, l in L))

    # === KKT conditions for leaders (x) ===
    @constraint(model,[i in 1:nl],
            [(- mu[1,i]*(0.2*y[1] - 1.333)/nl
            - mu[2,i]*(0.1*y[2] - 1)/nl
            - mu[3,i]*(-0.1)/nl
            - mu[4,i]*0.1/nl), x[i]] in MOI.Complements(2))

    # === Follower KKT conditions ===
    @constraint(model, [l in L],
        psi[1,l] == 2*(y[1]-3)
            - mu[1,l]*(1 + 0.2*xa + 2*y[4])
            - mu[3,l]*(0.333)
            - mu[4,l]*(-2*y[1])
            + xi[1,l]*s[1])

    @constraint(model, [l in L],
        psi[2,l] == 2*(y[2]-4)
            - mu[2,l]*(1 + 0.1*xa + 2*y[4])
            - mu[3,l]*(-1)
            - mu[4,l]*(-2*y[2])
            + xi[2,l]*s[2])

    @constraint(model, psi[3,1] == 0 - mu[1,1]*(-0.333) - mu[2,1] + xi[3,1]*s[3])
    @constraint(model, psi[3,2] == 2*(y[3]-1) - mu[1,2]*(-0.333) - mu[2,2] + xi[3,2]*s[3])
    @constraint(model, psi[3,3] == 0 - mu[1,3]*(-0.333) - mu[2,3] + xi[3,3]*s[3])
    @constraint(model, psi[3,4] == 2*(y[3]-1) - mu[1,4]*(-0.333) - mu[2,4] + xi[3,4]*s[3])

    @constraint(model, psi[4,1] == 0 - mu[1,1]*(2*y[1]) - mu[2,1]*(2*y[2]) + xi[4,1]*s[4])
    @constraint(model, psi[4,2] == 0 - mu[1,2]*(2*y[1]) - mu[2,2]*(2*y[2]) + xi[4,2]*s[4])
    @constraint(model, psi[4,3] == 20*y[4] - mu[1,3]*(2*y[1]) - mu[2,3]*(2*y[2]) + xi[4,3]*s[4])
    @constraint(model, psi[4,4] == 2*(y[4]-1) - mu[1,4]*(2*y[1]) - mu[2,4]*(2*y[2]) + xi[4,4]*s[4])

    @constraint(model, [i in I, l in L], sigma[i,l] == mu[i,l] + xi[i,l]*y[i])

    # === Followers' response constraints (with complementarity) ===
    @constraint(model, [l in L], s[1] == (1 + 0.2*xa)*y[1] - (3 + 1.333*xa) - 0.333*y[3] + 2*y[1]*y[4])
    @constraint(model, [l in L], s[2] == (1 + 0.1*xa)*y[2] - xa + y[3] + 2*y[2]*y[4])
    @constraint(model, [l in L], s[3] == 0.333*y[1] - y[2] + 1 - 0.1*xa)
    @constraint(model, [l in L], s[4] == 9 + 0.1*xa - y[1]^2 - y[2]^2)

    # === Complementarity constraints ===
    @constraint(model, [i in I, l in L], [y[i], psi[i,l]] in MOI.Complements(2))
    @constraint(model, [i in I, l in L], [s[i], sigma[i,l]] in MOI.Complements(2))
    @constraint(model, [i in I], [s[i], y[i]] in MOI.Complements(2))

    return model
end

#########################################################################
# outrata4-MPEC.mod
# Original AMPL coding by Sven Leyffer, Argonne National Laboratory, 2004.
#
# MPEC formulation of ...
#
# Multi-Leader-Follower Game (MLF)) derived from outrata31-outrata34,
# see J. Outrata, SIAM J. Optim. 4(2), pp.340ff, 1994. All Stackelberg
# players have the same constraints, but different objectives.
#
#########################################################################
function macepec_outrata4_model()
    nl = 4
    L = 1:nl
    I = 1:4

    model = Model()
    JuMP.set_name(model, "MacEPEC-outrata4")

    # ... Stackelberg leaders' variables
    @variable(model, 0 <= x[l in L] <= 10)
    @expression(model, xa, sum(x[l] for l in L) / nl)

    # ... Followers' variables
    @variable(model, y[i in I] >= 0)
    @variable(model, s[i in I] >= 0)

    # ... multipliers (differ for each leader)
    @variable(model, mF[l in L, i in I])
    @variable(model, my[l in L, i in I] >= 0)
    @variable(model, ms[l in L, i in I] >= 0)
    @variable(model, mC[l in L, i in I] >= 0)
    @variable(model, mx[l in L])

    # ... defined variables model rhs of KKT conditions
    @expression(model, kktx[l in L],
        mF[l,1]*(0.2*y[1] - 1.333) +
        mF[l,2]*(0.1*y[2] - 1) +
        mF[l,3]*(-0.1) +
        mF[l,4]*0.1 +
        mx[l]*nl
    )
    @expression(model, kkty1[l in L],
        mF[l,1]*(1 + 0.2*xa + 2*y[4]) +
        mF[l,3]*(0.333) +
        mF[l,4]*(-2*y[1]) +
        my[l,1] - mC[l,1]*s[1]
    )
    @expression(model, kkty2[l in L],
        mF[l,2]*(1 + 0.1*xa + 2*y[4]) +
        mF[l,3]*(-1) +
        mF[l,4]*(-2*y[2]) +
        my[l,2] - mC[l,2]*s[2]
    )
    @expression(model, kkty3[l in L],
        mF[l,1]*(-0.333) +
        mF[l,2] +
        my[l,3] - mC[l,3]*s[3]
    )
    @expression(model, kkty4[l in L],
        mF[l,1]*(2*y[1]) +
        mF[l,2]*(2*y[2]) +
        my[l,4] - mC[l,4]*s[4]
    )

    # === Objective ===
    @objective(model, Min, sum(mC[l, i] for l in L, i in I))

   # ... FO conditions of leader 1 (gradients in (x[1],y))
    @constraint(model, 2*(x[1] - 1) == kktx[1] / nl)
    @constraint(model, 2*(y[1] - 3) == kkty1[1])
    @constraint(model, 2*(y[2] - 4) == kkty2[1])
    @constraint(model, 2*(y[3] - 1) == kkty3[1])
    @constraint(model, 20*y[4] == kkty4[1])

    # ... FO conditions of leader 2 (gradients in (x[2],y))
    @constraint(model, 3*x[2] == kktx[2] / nl)
    @constraint(model, 2*(y[1] - 3) == kkty1[2])
    @constraint(model, 2*(y[2] - 4) == kkty2[2])
    @constraint(model, 2*(y[3] - 1) == kkty3[2])
    @constraint(model, 20*y[4] == kkty4[2])

    # ... FO conditions of leader 3 (gradients in (x[3],y))
    @constraint(model, 2*x[3] == kktx[3] / nl)
    @constraint(model, 2*(y[1] - 3) == kkty1[3])
    @constraint(model, 2*(y[2] - 4) == kkty2[3])
    @constraint(model, 2*(y[3] - 1) == kkty3[3])
    @constraint(model, 20*y[4] == kkty4[3])

    # ... FO conditions of leader 4 (gradients in (x[4],y))
    @constraint(model, 2*x[4] == kktx[4] / nl)
    @constraint(model, 2*(y[1] - 3) == kkty1[4])
    @constraint(model, 2*(y[2] - 4) == kkty2[4])
    @constraint(model, 2*(y[3] - 1) == kkty3[4])
    @constraint(model, 20*y[4] == kkty4[4])

    # ... first order conditions for the slacks (all players)
    @constraint(model, [l in L, i in I], 0 == -mF[l, i] + ms[l, i] - mC[l, i]*y[i])

    # ... followers reponse to average of leaders (xa)
    @constraint(model, [l in L],
        [s[1] - ((1 + 0.2*xa)*y[1] - (3 + 1.333*xa) - 0.333*y[3] + 2*y[1]*y[4]),
        mF[l,1]] in MOI.Complements(2))

    @constraint(model, [l in L],
        [s[2] - ((1 + 0.1*xa)*y[2] - xa + y[3] + 2*y[2]*y[4]),
        mF[l,2]] in MOI.Complements(2))

    @constraint(model, [l in L],
        [s[3] - (0.333*y[1] - y[2] + 1 - 0.1*xa),
        mF[l,3]] in MOI.Complements(2))

    @constraint(model, [l in L],
        [s[4] - (9 + 0.1*xa - y[1]^2 - y[2]^2),
        mF[l,4]] in MOI.Complements(2))

    # ... complementarity constraints
    @constraint(model, [l in L],
        [mx[l], x[l]] in MOI.Complements(2))
    @constraint(model, [l in L, i in I],
        [y[i], my[l, i]] in MOI.Complements(2))
    @constraint(model, [l in L, i in I],
        [s[i], ms[l, i]] in MOI.Complements(2))
    @constraint(model, [i in I], 0 >= s[i] * y[i])

    return model
end

