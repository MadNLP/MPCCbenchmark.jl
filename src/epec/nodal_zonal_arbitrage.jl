################################################################################
# Copyright (c) 2026: Quentin Lété
################################################################################

"""
    epec_nodal_zonal_arbitrage_model

Model implementing a nodal-zonal arbitrage problem using complementarity constraints.
The model comes from [1]. The model is applied to a simple 2-node 1-zone example sourced from [2].

## References
[1] "Power Generation Investment Under Zonal Electricity Pricing with Market-Based Re-dispatch", Lété et al. (2025)
[2] "Market-Based Redispatch in Zonal Electricity Markets: Inc-Dec Gaming as a Consequence of Inconsistent Power Market Design", Hirth & Slecht (2019)
"""
function epec_nodal_zonal_arbitrage_model()
    ## DATA DEFINITION
    # Generator capacities (100 units of 1GW each)
    X = ones(70)
    # Marginal costs for each generator
    MC = vcat(
        ones(20),          # 20GW of wind at 1€/MWh
        21.0:40.0,         # 20GW of COAL with increasing marginal cost
        66.0:70.0,         # 5GW of diesel peakers in the North
        41.0:65.0,
    )         # 30GW of diesel peakers in the South
    # Node of each generator: gas in node 2 (south), others in node 1 (north)
    Ng = vcat(ones(Int, 45), 2*ones(Int, 25))
    # Zone of each generator (all in zone 1)
    Zg = vcat(ones(Int, 70))
    # Generators per node
    Gn = [collect(1:45), collect(46:70)]
    # Generators per zone
    Gz = [collect(1:70)]
    # Nodes per zone
    Nz = [[1, 2]]
    # Zone per node
    Zn = [1, 1]
    # Transmission capacity (30GW)
    TC = [30]
    # Power Transfer Distribution Factor (PTDF)
    PTDF = [1 0;]
    # Nodal demand
    D = [0, 50]

    m = Model()

    ## VARIABLES
    @variable(m, y[i in 1:length(X)] >= 0)      # Generation in the day-ahead market
    @variable(m, y_tilde[1:length(X)])          # Redispatch generation
    @variable(m, mu[1:length(X)] >= 0)          # Dual variable for generation capacity
    @variable(m, mu_tilde[1:length(X)] >= 0)    # Dual variable for generation capacity
    @variable(m, delta[1:length(X)] >= 0)       # Dual variable for non-negativity of generation
    @variable(m, p[1:length(Gz)])               # Zonal power injection
    @variable(m, gamma)                         # Zonal price
    @variable(m, nu[1:length(Gn)])              # Nodal price adjustment
    @variable(m, r_tilde[1:length(Gn)])         # Redispatch nodal injection
    @variable(m, r[1:length(Gn)])               # Nodal power injection
    @variable(m, f[1:length(TC)])               # Flow on transmission lines
    @variable(m, psi[1:length(TC)])             # Dual variable for flow constraints
    @variable(m, phi)                           # Dual variable for power balance
    @variable(m, lambdap[1:length(TC)] >= 0)    # Dual variable for upper bound on flow
    @variable(m, lambdam[1:length(TC)] >= 0)    # Dual variable for lower bound on flow
    @variable(m, rho[1:length(Gz)])             # Zonal price
    @variable(m, rho_tilde[1:length(Gn)])       # Nodal price in redispatch

    ## COMPLEMENTARITY CONSTRAINTS
    @constraint(
        m,
        [g=1:length(X)],
        MC[g] + mu[g] + mu_tilde[g] - rho[Zg[g]] - delta[g] ⟂ y[g]
    ) # Generation complementarity
    @constraint(
        m,
        [g=1:length(X)],
        MC[g] + mu_tilde[g] - rho_tilde[Ng[g]] - delta[g] ⟂ y_tilde[g]
    ) # Redispatch generation complementarity
    @constraint(m, [g=1:length(X)], X[g] - y[g] - y_tilde[g] ⟂ mu_tilde[g]) # Redispatch capacity complementarity
    @constraint(m, [g=1:length(X)], X[g] - y[g] ⟂ mu[g]) # Day-ahead capacity complementarity
    @constraint(m, [g=1:length(X)], y[g] + y_tilde[g] ⟂ delta[g]) # Non-negative generation complementarity
    @constraint(m, [z=1:length(Gz)], rho[z] - gamma ⟂ p[z]) # Zonal power injection complementarity
    @constraint(m, 0 + sum(p[z] for z = 1:length(Gz)) ⟂ gamma) # Zonal balance constraint
    @constraint(
        m,
        [n=1:length(Gn)],
        r[n] - sum(y[g] for g in Gn[n]) + D[n] - r_tilde[n] ⟂ nu[n]
    ) # Nodal power balance complementarity
    @constraint(m, [n=1:length(Gn)], rho_tilde[n] + nu[n] ⟂ r_tilde[n]) # Redispatch nodal injection complementarity
    @constraint(
        m,
        [n=1:length(Gn)],
        -nu[n] - phi + sum(PTDF[l, n]*psi[l] for l = 1:length(TC)) ⟂ r[n]
    ) # Nodal power injection complementarity
    @constraint(m, [l=1:length(TC)], -psi[l] + lambdap[l] - lambdam[l] ⟂ f[l]) # Flow variable complementarity
    @constraint(
        m,
        [l=1:length(TC)],
        f[l] - sum(PTDF[l, n]*r[n] for n = 1:length(Gn)) ⟂ psi[l]
    ) # Flow constraint complementarity
    @constraint(m, sum(r[n] for n = 1:length(Gn)) ⟂ phi) # Power balance complementarity
    @constraint(m, [l=1:length(TC)], TC[l] - f[l] ⟂ lambdap[l]) # Upward transmission capacity complementarity
    @constraint(m, [l=1:length(TC)], TC[l] + f[l] ⟂ lambdam[l]) # Downward transmission capacity complementarity
    @constraint(
        m,
        [z=1:length(Gz)],
        -p[z] + sum(y[g] for g in Gz[z]) - sum(D[n] for n in Nz[z]) ⟂ rho[z]
    ) # Zonal power balance complementarity
    @constraint(
        m,
        [n=1:length(Gn)],
        -r_tilde[n] + sum(y_tilde[g] for g in Gn[n]) ⟂ rho_tilde[n]
    ) # Redispatch nodal power balance complementarity

    return m
end

function epec_nodal_zonal_arbitrage_large_model()
    data = JSON.parsefile(joinpath(@__DIR__, "data", "CWE_reduced.json"))

    # Unpack data
    G_Nidx = data["G_Nidx"]
    Gn = data["Gn"]
    G_Zidx = data["G_Zidx"]
    Gz = data["Gz"]
    DT = data["DT"]
    MC = data["MC"]
    IC = data["I"]
    FC = data["FC"]
    MC_ex = data["MC_ex"]
    IC_ex = data["IC_ex"]
    FC_ex = data["FC_ex"]
    VOLL = data["VOLL"]
    D = data["D"]
    N_names = data["N"]
    N = length(N_names)
    Z_names = data["Z"]
    Z = length(Z_names)
    Nz = data["Nz"]
    Zn = data["Zn"]
    X = data["X"]
    X_bar = data["X_bar"]
    L_names = data["L"]
    L = length(L_names)
    TC = data["TC"]
    PTDF = data["PTDF"]

    T = length(DT) # number of time periods
    NG = length(IC) # number of investable generators
    NG_ex = length(X_bar) # number of existing generators

    m = Model()
    # Investment
    @variable(m, x[i=1:NG, n=1:N] >= 0)          # new investments
    @variable(m, x_bar[g=1:NG_ex] >= 0)          # existing capacity chosen to keep
    # DA generation
    @variable(m, y[i=1:NG, t=1:T, n=1:N] >= 0)
    @variable(m, y_bar[g=1:NG_ex, t=1:T] >= 0)
    @variable(m, s[t=1:T, n=1:N] >= 0)           # load shedding
    # nodal injections
    @variable(m, r[t=1:T, n=1:N])
    @variable(m, p[t=1:T, z=1:Z])                # zonal injections
    @variable(m, pn[t=1:T, n=1:N])               # zonal injections
    # flows
    @variable(m, f[t=1:T, l=1:L])
    # redispatch
    @variable(m, y_tilde[i=1:NG, t=1:T, n=1:N])
    @variable(m, y_bar_tilde[g=1:NG_ex, t=1:T])
    @variable(m, s_tilde[t=1:T, n=1:N])
    @variable(m, r_tilde[t=1:T, n=1:N])
    @variable(m, f_tilde[t=1:T, l=1:L])
    # Prices
    @variable(m, rho[t=1:T, z=1:Z])              # zonal DA price
    @variable(m, rho_tilde[t=1:T, n=1:N])        # redispatch nodal price
    @variable(m, nu[t=1:T, n=1:N])               # nodal redispatch price
    # Network dual variables
    @variable(m, zeta[t=1:T, z=1:Z])
    @variable(m, phi[t=1:T])
    @variable(m, psi[t=1:T, l=1:L])
    @variable(m, lambdap[t=1:T, l=1:L] >= 0)
    @variable(m, lambdam[t=1:T, l=1:L] >= 0)
    @variable(m, phi_tilde[t=1:T])
    @variable(m, psi_tilde[t=1:T, l=1:L])
    @variable(m, lambdap_tilde[t=1:T, l=1:L] >= 0)
    @variable(m, lambdam_tilde[t=1:T, l=1:L] >= 0)
    # Duals for DA capacity
    @variable(m, mu[i=1:NG, t=1:T, n=1:N] >= 0)
    @variable(m, mu_bar[g=1:NG_ex, t=1:T] >= 0)
    @variable(m, mus[t=1:T, n=1:N] >= 0)
    @variable(m, mu_pos[i=1:NG, t=1:T, n=1:N] >= 0)
    @variable(m, mu_bar_pos[g=1:NG_ex, t=1:T] >= 0)
    @variable(m, mus_pos[t=1:T, n=1:N] >= 0)
    @variable(m, delta[i=1:NG] >= 0)
    @variable(m, delta_bar[g=1:NG_ex] >= 0)
    # Duals for redispatch
    @variable(m, mu_tilde[i=1:NG, t=1:T, n=1:N] >= 0)
    @variable(m, mus_tilde[t=1:T, n=1:N] >= 0)
    @variable(m, mu_bar_tilde[g=1:NG_ex, t=1:T] >= 0)

    # ---------- Investment -------------
    @constraint(
        m,
        [i=1:NG, n=1:N],
        IC[i] + FC[i] + delta[i] - sum(mu[i, t, n] + mu_tilde[i, t, n] for t = 1:T) ⟂
        x[i, n]
    )
    @constraint(
        m,
        [g=1:NG_ex],
        IC_ex[g] + FC_ex[g] + delta_bar[g] -
        sum(mu_bar[g, t] + mu_bar_tilde[g, t] for t = 1:T) ⟂ x_bar[g]
    )
    # ---------- DA generation ----------
    @constraint(
        m,
        [i=1:NG, t=1:T, n=1:N],
        DT[t]*MC[i] + mu[i, t, n] + mu_tilde[i, t, n] - rho[t, Zn[n]] - mu_pos[i, t, n] ⟂
        y[i, t, n]
    )
    @constraint(
        m,
        [g=1:NG_ex, t=1:T],
        DT[t]*MC_ex[g] + mu_bar[g, t] + mu_bar_tilde[g, t] - rho[t, G_Zidx[g]] -
        mu_bar_pos[g, t] ⟂ y_bar[g, t]
    )
    # ---------- redispatch ----------
    @constraint(
        m,
        [i=1:NG, t=1:T, n=1:N],
        DT[t]*MC[i] + mu_tilde[i, t, n] - rho_tilde[t, n] - mu_pos[i, t, n] ⟂
        y_tilde[i, t, n]
    )
    @constraint(
        m,
        [g=1:NG_ex, t=1:T],
        DT[t]*MC_ex[g] + mu_bar_tilde[g, t] - rho_tilde[t, G_Nidx[g]] - mu_bar_pos[g, t] ⟂
        y_bar_tilde[g, t]
    )
    # ---------- Capacity constraints ----------
    @constraint(m, [i=1:NG, t=1:T, n=1:N], x[i, n] - y[i, t, n] ⟂ mu[i, t, n])
    @constraint(
        m,
        [i=1:NG, t=1:T, n=1:N],
        x[i, n] - y[i, t, n] - y_tilde[i, t, n] ⟂ mu_tilde[i, t, n]
    )
    @constraint(m, [g=1:NG_ex, t=1:T], x_bar[g] - y_bar[g, t] ⟂ mu_bar[g, t])
    @constraint(
        m,
        [g=1:NG_ex, t=1:T],
        x_bar[g] - y_bar[g, t] - y_bar_tilde[g, t] ⟂ mu_bar_tilde[g, t]
    )
    @constraint(m, [t=1:T, n=1:N], D[n][t] - s[t, n] ⟂ mus[t, n])
    @constraint(m, [t=1:T, n=1:N], D[n][t] - s[t, n] - s_tilde[t, n] ⟂ mus_tilde[t, n])
    @constraint(m, [i=1:NG], X[i] - sum(x[i, n] for n = 1:N) ⟂ delta[i])
    @constraint(m, [g=1:NG_ex], X_bar[g] - x_bar[g] ⟂ delta_bar[g])
    # ---------- Non-negativity ----------
    @constraint(m, [i=1:NG, t=1:T, n=1:N], y[i, t, n] + y_tilde[i, t, n] ⟂ mu_pos[i, t, n])
    @constraint(m, [g=1:NG_ex, t=1:T], y_bar[g, t] + y_bar_tilde[g, t] ⟂ mu_bar_pos[g, t])
    @constraint(m, [t=1:T, n=1:N], s[t, n] + s_tilde[t, n] ⟂ mus_pos[t, n])
    # ---------- Load shedding ----------
    @constraint(
        m,
        [t=1:T, n=1:N],
        VOLL + mus[t, n] + mus_tilde[t, n] - rho[t, Zn[n]] ⟂ s[t, n]
    )
    @constraint(m, [t=1:T, n=1:N], VOLL + mus_tilde[t, n] - rho_tilde[t, n] ⟂ s_tilde[t, n])
    # ---------- Zonal balance ----------
    @constraint(m, [t=1:T, z=1:Z], rho[t, z] - zeta[t, z] ⟂ p[t, z])
    @constraint(m, [t=1:T, z=1:Z], p[t, z] - sum(pn[t, n] for n in Nz[z]) ⟂ zeta[t, z])
    # ---------- Nodal market clearing ----------
    @constraint(
        m,
        [t=1:T, n=1:N],
        r[t, n] - sum(y[i, t, n] for i = 1:NG) - sum(y_bar[g, t] for g in Gn[string(n)]) -
        s[t, n] + D[n][t] - r_tilde[t, n] ⟂ nu[t, n]
    )
    @constraint(m, [t=1:T, n=1:N], rho_tilde[t, n] + nu[t, n] ⟂ r_tilde[t, n])
    @constraint(
        m,
        [t=1:T, n=1:N],
        -nu[t, n] - phi_tilde[t] + sum(PTDF[l][n]*psi_tilde[t, l] for l = 1:L) ⟂ r[t, n]
    )
    @constraint(
        m,
        [t=1:T, n=1:N],
        zeta[t, Zn[n]] - phi[t] + sum(PTDF[l][n]*psi[t, l] for l = 1:L) ⟂ pn[t, n]
    )
    # ---------- PTDF & flow ----------
    @constraint(m, [t=1:T, l=1:L], -psi[t, l] + lambdap[t, l] - lambdam[t, l] ⟂ f[t, l])
    @constraint(
        m,
        [t=1:T, l=1:L],
        f[t, l] - sum(PTDF[l][n] * pn[t, n] for n = 1:N) ⟂ psi[t, l]
    )
    @constraint(m, [t=1:T], sum(pn[t, n] for n = 1:N) ⟂ phi[t])
    @constraint(
        m,
        [t=1:T, l=1:L],
        -psi_tilde[t, l] + lambdap_tilde[t, l] - lambdam_tilde[t, l] ⟂ f_tilde[t, l]
    )
    @constraint(
        m,
        [t=1:T, l=1:L],
        f_tilde[t, l] - sum(PTDF[l][n] * r[t, n] for n = 1:N) ⟂ psi_tilde[t, l]
    )
    @constraint(m, [t=1:T], sum(r[t, n] for n = 1:N) ⟂ phi_tilde[t])
    # ---------- Transmission bounds ----------
    @constraint(m, [t=1:T, l=1:L], TC[l] - f_tilde[t, l] ⟂ lambdap_tilde[t, l])
    @constraint(m, [t=1:T, l=1:L], TC[l] - f[t, l] ⟂ lambdap[t, l])
    @constraint(m, [t=1:T, l=1:L], TC[l] + f_tilde[t, l] ⟂ lambdam_tilde[t, l])
    @constraint(m, [t=1:T, l=1:L], TC[l] + f[t, l] ⟂ lambdam[t, l])
    # ---------- Zone to node redispatch ----------
    @constraint(
        m,
        [t=1:T, z=1:Z],
        -p[t, z] +
        sum(y[i, t, n] for i = 1:NG, n in Nz[z]) +
        sum(y_bar[g, t] for g in Gz[string(z)]) +
        sum(s[t, n] for n in Nz[z]) - sum(D[n][t] for n in Nz[z]) ⟂ rho[t, z]
    )
    @constraint(
        m,
        [t=1:T, n=1:N],
        -r_tilde[t, n] +
        sum(y_tilde[i, t, n] for i = 1:NG) +
        sum(y_bar_tilde[g, t] for g in Gn[string(n)]) +
        s_tilde[t, n] ⟂ rho_tilde[t, n]
    )

    return m
end
