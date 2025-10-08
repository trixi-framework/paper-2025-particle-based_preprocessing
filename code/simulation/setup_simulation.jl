using OrdinaryDiffEq
include(joinpath("..", "auxiliary_functions.jl"))

# ==========================================================================================
# ==== Experiment Setup
tspan = (0.0, 4.0)

particle_configuration = "packed"
free_slip_airfoil = true

# Boundary geometry and initial fluid particle positions
d_airfoil = 0.13895
chord_length = 1.0
domain_size = (4 * chord_length, 8 * d_airfoil)

factor = 25
particle_spacing = 1 / factor

const prescribed_velocity = 2.0
reynolds_number = 200
flow_direction = [1.0, 0.0]
fluid_density = 1.0
sound_speed = 20 * prescribed_velocity

material = (density=10 * fluid_density, E=1.4e3 * (fluid_density * prescribed_velocity^2),
            nu=0.4)

input_directory = joinpath(OUT_DIR, "airfoil", "NACA_6412")
geometry = load_geometry(joinpath(DATA_DIR, "airfoil", "NACA_6412.asc"))
solid_volume = TrixiParticles.volume(geometry)
fluid_volume = prod(domain_size) - solid_volume

fluid = vtk2trixi(joinpath(input_directory,
                           particle_configuration *
                           "_domain_dp_$(factor)_initial_condition.vtu"))
fluid.density .= fluid_density
fluid.mass .= fluid_density * fluid_volume / nparticles(fluid)
fluid.velocity[1, :] .= prescribed_velocity

solid_ = vtk2trixi(joinpath(input_directory,
                            particle_configuration *
                            "_dp_$(factor)_tlsph_initial_condition.vtu"))
solid_.density .= material.density
solid_.mass .= material.density * solid_volume / nparticles(solid_)

fixed_point = [0.04; 0.01;;]
fixed_point_radius = 0.025

# Find the fixed particles
if factor == 25
    candidates = particle_configuration == "lattice" ? [1, 2] : [29, 30]
else
    candidates = TrixiParticles.find_too_close_particles(solid_.coordinates, fixed_point,
                                                         fixed_point_radius)

    if factor == 50 && particle_configuration == "packed"
        push!(candidates, 22)
        push!(candidates, 67)
    end
end

fixed_particles = InitialCondition(; coordinates=solid_.coordinates[:, candidates],
                                   density=material.density,
                                   particle_spacing=particle_spacing)
ic_ = setdiff(solid_, fixed_particles)
solid = union(ic_, fixed_particles)

# open boundaries
pipe = RectangularTank(particle_spacing, domain_size, domain_size,
                       n_layers=4, faces=(false, false, true, true),
                       min_coordinates=(-1, -4 * d_airfoil), fluid_density)

open_boundary_layers = 6
n_buffer_particles = 5 * pipe.n_particles_per_dimension[2]^2
open_boundary_size = (particle_spacing * open_boundary_layers, domain_size[2])

min_coords_inlet = (-open_boundary_layers * particle_spacing, 0.0) .+ (-1, -4 * d_airfoil)
inlet = RectangularTank(particle_spacing, open_boundary_size, open_boundary_size,
                        fluid_density, n_layers=4, min_coordinates=min_coords_inlet,
                        faces=(false, false, true, true))

min_coords_outlet = (pipe.fluid_size[1], 0.0) .+ (-1, -4 * d_airfoil)
outlet = RectangularTank(particle_spacing, open_boundary_size, open_boundary_size,
                         fluid_density, n_layers=4, min_coordinates=min_coords_outlet,
                         faces=(false, false, true, true))

# ==========================================================================================
# ==== Fluid
smoothing_length = 1.3 * particle_spacing
smoothing_kernel = WendlandC2Kernel{2}()

fluid_density_calculator = ContinuityDensity()

kinematic_viscosity = prescribed_velocity * chord_length / reynolds_number

viscosity = ViscosityAdami(nu=kinematic_viscosity)

state_equation = StateEquationCole(; sound_speed, reference_density=fluid_density,
                                   exponent=1)

density_diffusion = DensityDiffusionMolteniColagrossi(delta=0.1)
fluid_system = WeaklyCompressibleSPHSystem(fluid, fluid_density_calculator,
                                           state_equation, smoothing_kernel,
                                           smoothing_length, viscosity=viscosity,
                                           buffer_size=n_buffer_particles,
                                           shifting_technique=TransportVelocityAdami(background_pressure=10_000),
                                           density_diffusion=density_diffusion)

# ==========================================================================================
# ==== Solid
solid_smoothing_length = smoothing_length # sqrt(2) * particle_spacing
solid_smoothing_kernel = smoothing_kernel # WendlandC2Kernel{2}()

# For the FSI we need the hydrodynamic masses and densities in the solid boundary model
hydrodynamic_densites = fluid_density * ones(size(solid.density))
hydrodynamic_masses = hydrodynamic_densites * particle_spacing^2

boundary_model_solid = BoundaryModelDummyParticles(hydrodynamic_densites,
                                                   hydrodynamic_masses,
                                                   state_equation=state_equation,
                                                   viscosity=free_slip_airfoil ? nothing :
                                                             viscosity,
                                                   AdamiPressureExtrapolation(),
                                                   smoothing_kernel, smoothing_length)

alpha_viscosity_tlsph = 0.05
viscosity_tlsph = ArtificialViscosityMonaghan(alpha=alpha_viscosity_tlsph)
solid_system = TotalLagrangianSPHSystem(solid,
                                        solid_smoothing_kernel, solid_smoothing_length,
                                        material.E, material.nu,
                                        n_clamped_particles=nparticles(fixed_particles),
                                        viscosity=viscosity_tlsph,
                                        penalty_force=PenaltyForceGanzenmueller(alpha=0.05),
                                        boundary_model=boundary_model_solid)

# ==========================================================================================
# ==== Open boundary
open_boundary_model = BoundaryModelMirroringTafuni(; mirror_method=ZerothOrderMirroring())

face_in = ([0.0, 0.0] .+ (-1, -4 * d_airfoil),
           [0.0, domain_size[2]] .+ (-1, -4 * d_airfoil))
inflow = BoundaryZone(; boundary_face=face_in, face_normal=flow_direction,
                      open_boundary_layers, density=fluid_density, particle_spacing,
                      boundary_type=InFlow(), initial_condition=inlet.fluid,
                      reference_velocity=(x, t) -> SVector(prescribed_velocity, 0.0))

face_out = ([domain_size[1], 0.0] .+ (-1, -4 * d_airfoil),
            [domain_size[1], domain_size[2]] .+ (-1, -4 * d_airfoil))
outflow = BoundaryZone(; boundary_face=face_out, face_normal=(-flow_direction),
                       open_boundary_layers, density=fluid_density, particle_spacing,
                       boundary_type=OutFlow(), initial_condition=outlet.fluid,
                       reference_pressure=0.0)

open_boundary = OpenBoundarySystem(inflow, outflow; fluid_system,
                                   boundary_model=open_boundary_model,
                                   buffer_size=n_buffer_particles)

# ==========================================================================================
# ==== Boundary

wall = union(pipe.boundary, inlet.boundary, outlet.boundary)

boundary_density_calculator = AdamiPressureExtrapolation()
viscosity_wall = nothing # no slip walls

boundary_model = BoundaryModelDummyParticles(wall.density, wall.mass,
                                             state_equation=state_equation,
                                             boundary_density_calculator,
                                             smoothing_kernel, smoothing_length,
                                             viscosity=viscosity_wall)

boundary_system = WallBoundarySystem(wall, boundary_model)

# ==========================================================================================
# ==== Simulation
min_corner = minimum(inlet.boundary.coordinates .- particle_spacing, dims=2)
max_corner = maximum(outlet.boundary.coordinates .+ particle_spacing, dims=2)

nhs = GridNeighborhoodSearch{2}(; cell_list=FullGridCellList(; min_corner, max_corner),
                                update_strategy=ParallelUpdate())

semi = Semidiscretization(fluid_system, solid_system, open_boundary,
                          boundary_system; neighborhood_search=nhs,
                          parallelization_backend=PolyesterBackend())
ode = semidiscretize(semi, tspan)

drag_force(system, data, t) = nothing
function drag_force(system::TotalLagrangianSPHSystem, data, t)
    dv = data.acceleration
    f_d = sum(TrixiParticles.eachparticle(system)) do particle
        return TrixiParticles.current_velocity(dv, system, particle) * system.mass[particle]
    end

    return f_d[1]
end

lift_force(system, data, t) = nothing
function lift_force(system::TotalLagrangianSPHSystem, data, t)
    dv = data.acceleration
    f_l = sum(TrixiParticles.eachparticle(system)) do particle
        return TrixiParticles.current_velocity(dv, system, particle) * system.mass[particle]
    end

    return f_l[2]
end

marker_point = [0.85; 0.03;;]
deflection_candidates = TrixiParticles.find_too_close_particles(solid.coordinates,
                                                                marker_point,
                                                                1.3 * particle_spacing)

initial_center = Ref(SVector(zero(eltype(solid_system)), zero(eltype(solid_system))))
for particle in deflection_candidates
    initial_center[] += TrixiParticles.initial_coords(solid_system, particle)
end
initial_center[] /= length(deflection_candidates)

deflection_y(system, data, t) = nothing
function deflection_y(system::TotalLagrangianSPHSystem, data, t)
    current_center = SVector(zero(eltype(system)), zero(eltype(system)))

    for particle in deflection_candidates
        current_center += TrixiParticles.current_coords(system, particle)
    end

    current_center /= length(deflection_candidates)

    displacement = current_center - initial_center[]

    return displacement[2] + initial_center[][2]
end

output_directory = joinpath(OUT_DIR, "airfoil", "NACA_6412",
                            "dp_$(factor)_tvf" * (free_slip_airfoil ? "_free_slip" : ""))

info_callback = InfoCallback(interval=100)
saving_callback = SolutionSavingCallback(dt=0.02, prefix=particle_configuration,
                                         output_directory=output_directory)
extra_callback = nothing
pp_cb_ekin = PostprocessCallback(; dt=0.02, deflection_y=deflection_y,
                                 f_l=lift_force, f_d=drag_force,
                                 filename=particle_configuration * "_results",
                                 write_file_interval=1, output_directory=output_directory)

callbacks = CallbackSet(info_callback, saving_callback, UpdateCallback(), extra_callback,
                        pp_cb_ekin)

# Use a Runge-Kutta method with automatic (error based) time step size control.
# Limiting of the maximum stepsize is necessary to prevent crashing.
# When particles are approaching a wall in a uniform way, they can be advanced
# with large time steps. Close to the wall, the stepsize has to be reduced drastically.
# Sometimes, the method fails to do so because forces become extremely large when
# fluid particles are very close to boundary particles, and the time integration method
# interprets this as an instability.
sol = solve(ode, RDPK3SpFSAL35(),
            abstol=1e-8, # Default abstol is 1e-6 (may need to be tuned to prevent boundary penetration)
            reltol=1e-4, # Default reltol is 1e-3 (may need to be tuned to prevent boundary penetration)
            dtmax=1e-2, # Limit stepsize to prevent crashing
            save_everystep=false, callback=callbacks);
