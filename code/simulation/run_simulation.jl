include(joinpath("..", "auxiliary_functions.jl"))

free_slip = true

# factor_resolution = 25
trixi_include(joinpath(CODE_DIR, "simulation", "setup_simulation.jl"),
              free_slip_airfoil=free_slip, particle_configuration="lattice",
              alpha_viscosity_tlsph=0.05, factor=25)
trixi_include(joinpath(CODE_DIR, "simulation", "setup_simulation.jl"),
              free_slip_airfoil=free_slip, particle_configuration="packed",
              alpha_viscosity_tlsph=0.05, factor=25)

# factor_resolution = 50
trixi_include(joinpath(CODE_DIR, "simulation", "setup_simulation.jl"),
              free_slip_airfoil=free_slip, particle_configuration="lattice",
              alpha_viscosity_tlsph=0.05, factor=50)
trixi_include(joinpath(CODE_DIR, "simulation", "setup_simulation.jl"),
              free_slip_airfoil=free_slip, particle_configuration="packed",
              alpha_viscosity_tlsph=0.05, factor=50)

# factor_resolution = 100
trixi_include(joinpath(CODE_DIR, "simulation", "setup_simulation.jl"),
              free_slip_airfoil=free_slip, particle_configuration="lattice",
              alpha_viscosity_tlsph=0.1, factor=100)
trixi_include(joinpath(CODE_DIR, "simulation", "setup_simulation.jl"),
              free_slip_airfoil=free_slip, particle_configuration="packed",
              alpha_viscosity_tlsph=0.1, factor=100)

# factor_resolution = 200
trixi_include(joinpath(CODE_DIR, "simulation", "setup_simulation.jl"),
              free_slip_airfoil=free_slip, particle_configuration="lattice",
              alpha_viscosity_tlsph=0.2, factor=200)
trixi_include(joinpath(CODE_DIR, "simulation", "setup_simulation.jl"),
              free_slip_airfoil=free_slip, particle_configuration="packed",
              alpha_viscosity_tlsph=0.2, factor=200)

# factor_resolution = 400
trixi_include(joinpath(CODE_DIR, "simulation", "setup_simulation.jl"),
              free_slip_airfoil=free_slip, particle_configuration="lattice",
              alpha_viscosity_tlsph=0.4, factor=400)
trixi_include(joinpath(CODE_DIR, "simulation", "setup_simulation.jl"),
              free_slip_airfoil=free_slip, particle_configuration="packed",
              alpha_viscosity_tlsph=0.4, factor=400)
