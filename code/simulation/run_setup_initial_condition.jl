include(joinpath("..", "auxiliary_functions.jl"))

# factor_resolution = 25
trixi_include(joinpath(CODE_DIR, "simulation", "setup_initial_condition.jl"), factor=25)

# factor_resolution = 50
trixi_include(joinpath(CODE_DIR, "simulation", "setup_initial_condition.jl"), factor=50)

# factor_resolution = 100
trixi_include(joinpath(CODE_DIR, "simulation", "setup_initial_condition.jl"), factor=100)

# factor_resolution = 200
trixi_include(joinpath(CODE_DIR, "simulation", "setup_initial_condition.jl"), factor=200)

# factor_resolution = 400
trixi_include(joinpath(CODE_DIR, "simulation", "setup_initial_condition.jl"), factor=400)
