
function reformulate_to_vertical!(model::MOI.ModelLike)
    VF = Union{MOI.VectorAffineFunction,MOI.VectorQuadraticFunction,MOI.VectorNonlinearFunction}
    contypes = MOI.get(model, MOI.ListOfConstraintTypesPresent())
    for (F, S) in contypes
        # Parse only complementarity constraints
        if F <: VF && S <: MOI.Complements
            conindices = MOI.get(model, MOI.ListOfConstraintIndices{F,S}())
            for cidx in conindices
                fun = MOI.get(model, MOI.ConstraintFunction(), cidx)
                set = MOI.get(model, MOI.ConstraintSet(), cidx)
                ComplementOpt.reformulate_to_vertical!(model, fun, set)
                MOI.delete(model, cidx)
            end
        end
    end
    return model
end

function reformulate_to_nonlinear!(model, relaxation)
    cc_cons = MOI.get(
        model,
        MOI.ListOfConstraintIndices{MOI.VectorOfVariables, MOI.Complements}(),
    )
    for cidx in cc_cons
        fun = MOI.get(model, MOI.ConstraintFunction(), cidx)
        set = MOI.get(model, MOI.ConstraintSet(), cidx)
        ComplementOpt.reformulate_as_nonlinear_program!(model, relaxation, fun, set)
        MOI.delete(model, cidx)
    end
    return model
end

"""
    reformulate_to_standard_form!(model::ModelLike)

Reformulate the mixed-complementarity constraints ``x1 ⟂ (lb <= x2 <= ub)`` with slack variables.

If `support_lb=true`, the mixed-complementarity constraints are reformulated as:
```
x2n = ub - x2
x1 = x1p - x1n
0 <= x1p ⟂ x2  >= lb
0 <= x1n ⟂ x2n >= 0

```
Otherwise, if `support_lb=false`, a slack `x2p` is added to the problem:
```
x2p = x2 - lb
x2n = ub - x2
x1 = x1p - x1n
0 <= x1p ⟂ x2p >= 0
0 <= x1n ⟂ x2n >= 0

```

The function assumes the model is passed in vertical form, and returns an error otherwise.
"""
function reformulate_to_standard_form!(model::MOI.ModelLike; support_lb=true)
    cc_cons = MOI.get(model, MOI.ListOfConstraintIndices{MOI.VectorOfVariables, MOI.Complements}())
    ind_cc1, ind_cc2 = MOI.VariableIndex[], MOI.VariableIndex[]

    for cidx in cc_cons
        fun = MOI.get(model, MOI.ConstraintFunction(), cidx)
        set = MOI.get(model, MOI.ConstraintSet(), cidx)
        MOI.delete(model, cidx)
        n_comp = div(set.dimension, 2)
        for cc in 1:n_comp
            x1 = fun.variables[cc]
            x2 = fun.variables[cc + n_comp]
            # Get bounds on x1
            lb1, ub1 = MOIU.get_bounds(model, Float64, x1)
            # Get bounds on x2
            lb2, ub2 = MOIU.get_bounds(model, Float64, x2)

            # If lb2 = ub2, the left-hand-side x1 is free and we discard this complementarity constraint.
            if lb2 == ub2
                continue
            end

            if isfinite(lb2) && isinf(ub2)
                # Ensure x1 is well posed
                if !iszero(lb1) || isfinite(ub1)
                    error("The problem does not follow MOI convention for mixed-complementarity constraints: the LHS lower bound is not zero.")
                end
                if isinf(lb1)
                    MOI.add_constraint(model, x1, MOI.GreaterThan(0.0))
                end
                # Add a slack x2p = x2 - lb2 if lb2 ≠ 0
                if !support_lb && !iszero(lb2)
                    x2p, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
                    MOI.add_constraint(
                        model,
                        1.0 * x2 - 1.0 * x2p,
                        MOI.EqualTo(lb2),
                    )
                    push!(ind_cc1, x1)
                    push!(ind_cc2, x2p)
                else
                    # Nothing to change
                    push!(ind_cc1, x1)
                    push!(ind_cc2, x2)
                end
            elseif isinf(lb2) && isfinite(ub2)
                # Ensure x1 is well posed
                if !iszero(ub1) || isfinite(lb1)
                    error("The problem does not follow MOI convention for mixed-complementarity constraints: the LHS upper bound is not zero.")
                end
                # Add a slack x1n = -x1
                x1n, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
                # x1n = -x1
                MOI.add_constraint(
                    model,
                    1.0 * x1 + 1.0 * x1n,
                    MOI.EqualTo(0.0),
                )
                # Add a slack x2n = ub2 - x2
                x2n, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
                MOI.add_constraint(
                    model,
                    1.0*x2 + 1.0*x2n,
                    MOI.EqualTo(ub2),
                )
                push!(ind_cc1, x1n)
                push!(ind_cc2, x2n)
            elseif isfinite(lb2) && isfinite(ub2)
                # Reformulate the mixed-complementarity constraint as two complementarity constraints.
                # First, ensure x1 is well posed
                if isfinite(lb1) || isfinite(ub1)
                    error("The problem does not follow MOI convention for mixed-complementarity constraints: the LHS should be free.")
                end

                # Add two slacks such that x1 = x1p - x1n
                x1p, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
                x1n, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
                MOI.add_constraint(
                    model,
                    1.0*x1 - 1.0*x1p + 1.0*x1n,
                    MOI.EqualTo(0.0),
                )
                push!(ind_cc1, x1p)
                if !support_lb && !iszero(lb2)
                    # Add a slack x2p = x2 - lb2 if lb2 ≠ 0
                    x2p, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
                    MOI.add_constraint(
                        model,
                        1.0*x2 - 1.0*x2p,
                        MOI.EqualTo(lb2),
                    )
                    push!(ind_cc2, x2p)
                else
                    push!(ind_cc2, x2)
                end

                x2n, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
                MOI.add_constraint(
                    model,
                    1.0*x2 + 1.0*x2n,
                    MOI.EqualTo(ub2),
                )
                push!(ind_cc1, x1n)
                push!(ind_cc2, x2n)
            else
                error("Problem is not well specified")
            end
        end
    end
    @assert length(ind_cc1) == length(ind_cc2)
    return ind_cc1, ind_cc2
end

