
include("collocation.jl")

struct NOSNOCBenchmark <: AbstractBenchmarkSetting
    collocation::RKScheme
end

get_tag(::NOSNOCBenchmark) = "nosnoc"

include("acrobot_model.jl")
include("bilinear_spring_damper_model.jl")
include("cart_pole_with_friction.jl")
include("liquid_gas_tank_model.jl")
include("ocp_motor_with_friction.jl")
include("schumacher.jl")
include("sliding_mode_ocp.jl")
include("stewart_anitescu.jl")

function load_model(instance, bench::NOSNOCBenchmark)
    func = instance[1]
    return func(instance[2]..., bench.collocation)
end

function get_name(instance, ::NOSNOCBenchmark)
    func = instance[1]
    id = match(r"^nosnoc_(.*)_model", string(func)).captures[1]
    k = prod(instance[2])
    return "$(id)_$(k)"
end

function get_instances(bench::NOSNOCBenchmark)
    return [
        (nosnoc_acrobot_model, (30, 2)),
        (nosnoc_acrobot_model, (50, 2)),
        (nosnoc_acrobot_model, (60, 2)),
        (nosnoc_bilinear_spring_damper_model, (30, 2)),
        (nosnoc_bilinear_spring_damper_model, (44, 2)),
        (nosnoc_bilinear_spring_damper_model, (60, 2)),
        (nosnoc_cart_pole_with_friction_model, (30, 2)),
        (nosnoc_cart_pole_with_friction_model, (33, 2)),
        (nosnoc_cart_pole_with_friction_model, (50, 2)),
        (nosnoc_cart_pole_with_friction_model, (57, 2)),
        (nosnoc_liquid_gas_tank_model, (50, 2)),
        (nosnoc_liquid_gas_tank_model, (100, 2)),
        (nosnoc_motor_with_friction_model, (27, 3)),
        (nosnoc_motor_with_friction_model, (30, 3)),
        (nosnoc_motor_with_friction_model, (50, 3)),
        (nosnoc_motor_with_friction_model, (53, 3)),
        (nosnoc_schumacher_model, (50, 2)),
        (nosnoc_schumacher_model, (80, 2)),
        (nosnoc_schumacher_model, (120, 2)),
        (nosnoc_sliding_mode_ocp_model, (6, 3)),
        (nosnoc_sliding_mode_ocp_model, (37, 3)),
        (nosnoc_sliding_mode_ocp_model, (50, 3)),
    ]
end

