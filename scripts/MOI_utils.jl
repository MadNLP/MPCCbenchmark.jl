
reformulate_to_vertical!(model::JuMP.Model) = reformulate_to_vertical!(JuMP.backend(model))
function reformulate_to_vertical!(model::MOI.ModelLike)
    ind_cc1, ind_cc2 = MOI.VariableIndex[], MOI.VariableIndex[]

    contypes = MOI.get(model, MOI.ListOfConstraintTypesPresent())
    for (F, S) in contypes
        # Parse only complementarity constraints
        if S == MOI.Complements
            conindices = MOI.get(model, MOI.ListOfConstraintIndices{F, S}())
            for cidx in conindices
                fun = MOI.get(model, MOI.ConstraintFunction(), cidx)
                set = MOI.get(model, MOI.ConstraintSet(), cidx)
                n_comp = div(set.dimension, 2)
                if isa(fun, MOI.VectorOfVariables)
                    append!(ind_cc1, fun.variables[1:n_comp])
                    append!(ind_cc2, fun.variables[n_comp+1:end])
                elseif isa(fun, VF)
                    # Read each complementarity constraint and get corresponding indices
                    cc_lhs, cc_rhs = ComplementOpt._parse_complementarity_constraint(fun, n_comp)
                    for (lhs, x2) in zip(cc_lhs, cc_rhs)
                        # Check if x2 is bounded.
                        lb, ub = MOIU.get_bounds(model, Float64, x2)
                        if isinf(lb) && isinf(ub)
                            # If x2 is unbounded, the LHS is directly converted to an equality constraint.
                            MOI.add_constraint(model, lhs, MOI.EqualTo{Float64}(0))
                        elseif isa(lhs, MOI.VariableIndex)
                            # If lhs is a variable, no need to reformulate the
                            # complementarity constraint using a slack.
                            # TODO: we should check if the variable lhs is bounded.
                            push!(ind_cc1, lhs)
                            push!(ind_cc2, x2)
                        else
                            # Else, reformulate LHS using vertical form
                            x1 = MOI.add_variable(model)
                            new_lhs = MOIU.operate!(-, Float64, lhs, x1)
                            MOI.add_constraint(model, new_lhs, MOI.EqualTo{Float64}(0))
                            push!(ind_cc1, x1)
                            push!(ind_cc2, x2)
                        end
                    end
                else
                    error("Complementary constraints formulated with $(typeof(fun)) are not yet supported")
                end
                # We delete the complementarity constraints
                MOI.delete(model, cidx)
            end
        end
    end
    return ind_cc1, ind_cc2
end

