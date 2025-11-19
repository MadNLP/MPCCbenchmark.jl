
function compute_branch_limits!(data; default_pad=1.0472)
    nl = length(data.branch)
    for i in 1:length(data.branch)
        b = data.branch[i]

        if b.rate_a <= 0

            angmin = b.angmin
            if iszero(angmin) || b.angmin <= -pi/2
                angmin = -default_pad
            end
            angmax = b.angmax
            if  b.angmax >= pi/2
                angmax = default_pad
            end
            if iszero(angmin) && iszero(angmax)
                angmin = -default_pad
                angmax =  default_pad
            end

            theta_max = max(abs(angmin), abs(angmax))

            r = b.br_r
            x = b.br_x
            z = r + im * x
            y = LinearAlgebra.pinv(z)
            y_mag = abs.(y)

            fr_vmax = data.bus[b.f_bus].vmax
            to_vmax = data.bus[b.t_bus].vmax
            m_vmax = max(fr_vmax, to_vmax)

            c_max = sqrt(fr_vmax^2 + to_vmax^2 - 2*fr_vmax*to_vmax*cos(theta_max))
            rate_a = y_mag*m_vmax*c_max

            data.branch[i] = ExaPowerIO.BranchData{Float64}(
                b.i,
                b.f_bus,
                b.t_bus,
                b.br_r,
                b.br_x,
                b.b_fr,
                b.b_to,
                b.g_fr,
                b.g_to,
                rate_a,
                b.rate_b,
                b.rate_c,
                b.tap,
                b.shift,
                b.status,
                angmin,
                angmax,
                b.f_idx,
                b.t_idx
            )
            data.arc[i] = ExaPowerIO.ArcData(i, b.f_bus, rate_a)
            data.arc[i+nl] = ExaPowerIO.ArcData(i+nl, b.t_bus, rate_a)
        end
    end
end

