
function get_complementarity_residual(model::JuMP.Model, ind_cc1, ind_cc2)
    moi_model = JuMP.backend(model)
    x1 = MOI.get.(moi_model, MOI.VariablePrimal(), ind_cc1)
    x2 = MOI.get.(moi_model, MOI.VariablePrimal(), ind_cc2)

    bounds = MOI.Utilities.get_bounds.(moi_model, Float64, ind_cc2)
    lb = [b[1] for b in bounds]
    ub = [b[2] for b in bounds]

    resid = max.(min.(x1, x2 .- lb), x2 .- ub)
    return norm(resid, Inf)
end

function get_complementarity_residual(nlp::AbstractNLPModel, solution::AbstractVector, ind_cc1, ind_cc2)
    lb = NLPModels.get_lvar(nlp)
    ub = NLPModels.get_uvar(nlp)
    x1 = solution[ind_cc1]
    x2 = solution[ind_cc2]
    resid = min.(x1 .- lb[ind_cc1], x2 .- lb[ind_cc2])
    return norm(resid, Inf)
end

