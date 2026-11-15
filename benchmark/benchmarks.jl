using BenchmarkTools
using LinearSolve: QRFactorization
using Magnetostatics
using StaticArrays

function eval_vec_batch(f::F, points::Vector{SVector{3, Float64}}) where {F}
    acc = zero(SVector{3, Float64})
    @inbounds for p in points
        acc += f(p)
    end
    return acc
end

function eval_vec_batch(f::F, arg::A, points::Vector{SVector{3, Float64}}) where {F, A}
    acc = zero(SVector{3, Float64})
    @inbounds for p in points
        acc += f(arg, p)
    end
    return acc
end

function eval_scalar_batch(
        f::F, arg::A, points::Vector{SVector{3, Float64}}
    ) where {F, A}
    acc = 0.0
    @inbounds for p in points
        acc += f(arg, p)
    end
    return acc
end

function eval_solve_batch(
        solver::S, source::C, points::Vector{SVector{3, Float64}}
    ) where {S, C}
    acc = zero(SVector{3, Float64})
    @inbounds for p in points
        acc += solve(solver, source, p)
    end
    return acc
end

function eval_loop_batch(loop::L, points::Vector{SVector{3, Float64}}) where {L}
    acc = zero(SVector{3, Float64})
    @inbounds for p in points
        acc += getB_loop(p, loop)
    end
    return acc
end

const SUITE = BenchmarkGroup()

SUITE["solvers"] = BenchmarkGroup()

# Biot-Savart
SUITE["solvers"]["biot_savart"] = BenchmarkGroup()
solver, wire, biot_points = let
    points = vec(
        [
            SVector(x, y, z) for x in (-0.45, 0.45), y in (-0.45, 0.45),
                z in (-0.35, -0.1, 0.1, 0.35)
        ]
    )
    loop = CurrentLoop(1.0, 1.0, SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 1.0))
    # creates a 100-segment wire
    wire = discretize_loop(loop, 100)
    solver = BiotSavart()

    solver, wire, points
end

SUITE["solvers"]["biot_savart"]["discretized_loop"] =
    @benchmarkable eval_solve_batch($solver, $wire, $biot_points)

# FFT Solver
SUITE["solvers"]["fft"] = BenchmarkGroup()
fft_solver, J, dx = let
    solver = FFTSolver()
    dx = 0.1
    # 3x32x32x32 array
    J = rand(3, 32, 32, 32)
    solver, J, dx
end

SUITE["solvers"]["fft"]["small_grid"] = @benchmarkable solve($fft_solver, $J, $dx)
fft_open_solver = FFTSolver(OpenBoundary())
SUITE["solvers"]["fft"]["open_grid"] = @benchmarkable solve($fft_open_solver, $J, $dx)

# Vector Potential
SUITE["solvers"]["potential"] = BenchmarkGroup()
loop_points = vec(
    [
        SVector(x, y, z) for x in (-0.75, -0.35, -0.15, 0.15, 0.35),
            y in (-0.35, -0.15, 0.15, 0.35, 0.75), z in (-0.4, -0.15, 0.15, 0.4)
    ]
)
pot_solver, pot_loop = let
    solver = VectorPotential()
    loop = CurrentLoop(1.0, 1.0, SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 1.0))
    solver, loop
end

SUITE["solvers"]["potential"]["current_loop"] =
    @benchmarkable eval_solve_batch($pot_solver, $pot_loop, $loop_points)

# Analytical
SUITE["solvers"]["analytical"] = BenchmarkGroup()
loop = CurrentLoop(3.0, 1.0, SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 1.0))

SUITE["solvers"]["analytical"]["current_loop"] =
    @benchmarkable eval_loop_batch($loop, $loop_points)

dipole = Dipole(SVector(0.0, 0.0, 1.0))
dipole_points = vec(
    [
        SVector(x, y, z) for x in (-0.75, -0.25, 0.25, 0.75),
            y in range(-0.7, 0.7; length = 8), z in range(-0.7, 0.7; length = 8)
    ]
)

SUITE["solvers"]["analytical"]["dipole"] =
    @benchmarkable eval_vec_batch($dipole, $dipole_points)

solenoid_thin = Solenoid(0.5, 1.0, 1000.0)
solenoid_thick = Solenoid(0.5, 1.0, 1000.0; thickness = 0.1)
tokamak_coils = TokamakCoils(1.0, 1.0, 1.0e3, 1.0e3)
cuboid_mag = CuboidMagnet(SVector(0.2, 0.3, 0.4), SVector(0.0, 0.0, 1.0e6))
thick_points = vec(
    [
        SVector(x, y, z) for x in (-0.25, 0.25), y in (-0.25, 0.25),
            z in (-0.35, -0.1, 0.1, 0.35)
    ]
)
magnet_points = vec(
    [
        SVector(x, y, z) for x in (-0.35, -0.05, 0.05, 0.35),
            y in (-0.45, -0.05, 0.05, 0.45), z in (-0.55, -0.05, 0.05, 0.55)
    ]
)
SUITE["solvers"]["analytical"]["solenoid_thin"] =
    @benchmarkable eval_vec_batch($solenoid_thin, $loop_points)
SUITE["solvers"]["analytical"]["solenoid_thick"] =
    @benchmarkable eval_vec_batch($solenoid_thick, $thick_points)
SUITE["solvers"]["analytical"]["tokamak_coils"] =
    @benchmarkable eval_vec_batch($tokamak_coils, $thick_points)
SUITE["solvers"]["analytical"]["cuboid_magnet"] =
    @benchmarkable eval_vec_batch($cuboid_mag, $magnet_points)

# Field Composition
SUITE["solvers"]["composition"] = BenchmarkGroup()
composed_field = let
    rot = SMatrix{3, 3, Float64}(
        cos(0.3), 0.0, -sin(0.3),
        0.0, 1.0, 0.0,
        sin(0.3), 0.0, cos(0.3)
    )
    2.0 * translate_field(dipole, SVector(0.1, -0.2, 0.05)) + rotate_field(cuboid_mag, rot)
end
SUITE["solvers"]["composition"]["transformed_sum"] =
    @benchmarkable eval_vec_batch($composed_field, $magnet_points)

# Poisson Solver
SUITE["solvers"]["poisson"] = BenchmarkGroup()
poisson_solver, poisson_dst_solver, poisson_J = let
    xs = range(0.0, 1.0; length = 16)
    ys = range(0.0, 1.0; length = 16)
    zs = range(0.0, 1.0; length = 16)
    solver = PoissonSolver(xs, ys, zs)
    dst_solver = PoissonSolver(xs, ys, zs; alg = DSTFactorization())
    J = rand(3, 16, 16, 16)
    solver, dst_solver, J
end

SUITE["solvers"]["poisson"]["small_grid"] =
    @benchmarkable solve($poisson_solver, $poisson_J)
SUITE["solvers"]["poisson"]["dst_grid"] =
    @benchmarkable solve($poisson_dst_solver, $poisson_J)

# Interpolated Grid Magnetic Field
SUITE["solvers"]["grid_field"] = BenchmarkGroup()
grid_field, grid_points = let
    gs = range(0.0, 1.0; length = 16)
    gf = GridMagneticField(rand(3, 16, 16, 16), gs, gs, gs)
    points = vec(
        [
            SVector(x, y, z) for x in (0.15, 0.35, 0.65, 0.85),
                y in range(0.1, 0.9; length = 4), z in range(0.12, 0.88; length = 8)
        ]
    )
    gf, points
end
SUITE["solvers"]["grid_field"]["eval_B"] =
    @benchmarkable eval_vec_batch($grid_field, $grid_points)
SUITE["solvers"]["grid_field"]["eval_J"] =
    @benchmarkable eval_vec_batch(current_density, $grid_field, $grid_points)

# Axisymmetric field reconstruction from sparse probes: FLARE coils sampled on a
# 12 × 11 probe array, with virtual filaments for the coil field and a loose
# minimum-current prior on a collocation grid.
SUITE["reconstruction"] = BenchmarkGroup()
recon_problem = let
    flare = FLARECoils(; I_pf = 45.0e3, I_ef = 12.0e3, B_guide = 0.0)
    domain = (0.45, 1.25, -0.28, 0.28)
    probes = [
        (R, Z) for R in range(0.48, 1.22; length = 12)
            for Z in range(-0.25, 0.25; length = 11)
    ]
    positions = [p[i] for p in probes, i in 1:2]
    B = reduce(vcat, [flare(SVector(R, 0.0, Z))' for (R, Z) in probes])
    grid(h) = [
        p[i] for p in [
                (R, Z) for R in domain[1]:h:domain[2]
                for Z in domain[3]:h:domain[4]
            ], i in 1:2
    ]
    Rmin, Rmax, Zmin, Zmax = domain .+ (-0.1, 0.1, -0.1, 0.1)
    perimeter = 2 * (Rmax - Rmin + Zmax - Zmin)
    filaments = map(range(0, perimeter; length = 49)[1:(end - 1)]) do s
        s < Rmax - Rmin && return (Rmin + s, Zmin)
        s -= Rmax - Rmin
        s < Zmax - Zmin && return (Rmax, Zmin + s)
        s -= Zmax - Zmin
        s < Rmax - Rmin && return (Rmax - s, Zmax)
        return (Rmin, Zmax - (s - (Rmax - Rmin)))
    end
    collocation = grid(0.015)
    common = (;
        domain, regularization = 1.0e2,
        uncertainties = fill(5.0e-4, size(B)),
        filament_positions = [f[i] for f in filaments, i in 1:2],
        current_positions = collocation,
        current_density = zeros(size(collocation, 1)),
        current_uncertainties = fill(1.0e5, size(collocation, 1)),
    )
    (; positions, B, grid, common)
end

SUITE["reconstruction"]["fit"] = BenchmarkGroup()
SUITE["reconstruction"]["fit"]["wendland_lsmr"] = @benchmarkable(
    reconstruct_axisymmetric_field(
        p.positions, p.B; length_scale = 0.12, basis_centers = p.grid(0.06),
        p.common...
    ),
    setup = (p = $recon_problem)
)
SUITE["reconstruction"]["fit"]["gaussian_qr"] = @benchmarkable(
    reconstruct_axisymmetric_field(
        p.positions, p.B; kernel = :gaussian, length_scale = 0.06,
        basis_centers = p.grid(0.03), algorithm = QRFactorization(), p.common...
    ),
    setup = (p = $recon_problem)
)
if hasfield(AxisymmetricFieldReconstruction, :element_currents)
    SUITE["reconstruction"]["fit"]["wendland_elements_qr"] = @benchmarkable(
        reconstruct_axisymmetric_field(
            p.positions, p.B; length_scale = 0.12, basis_centers = p.grid(0.06),
            current_element_positions = p.grid(0.05),
            algorithm = QRFactorization(), p.common...
        ),
        setup = (p = $recon_problem)
    )
end

recon_field = reconstruct_axisymmetric_field(
    recon_problem.positions, recon_problem.B; length_scale = 0.12,
    basis_centers = recon_problem.grid(0.06), recon_problem.common...
)
recon_points = vec(
    [
        SVector(R * cos(phi), R * sin(phi), z) for R in range(0.52, 1.18; length = 4),
            phi in (0.1, 0.4), z in range(-0.22, 0.22; length = 4)
    ]
)
recon_j_points = vec(
    [
        SVector(R * cos(phi), R * sin(phi), z) for R in range(0.52, 1.18; length = 4),
            phi in (0.05, 0.25, 0.45, 0.65), z in range(-0.22, 0.22; length = 4)
    ]
)

SUITE["reconstruction"]["evaluate"] = BenchmarkGroup()
SUITE["reconstruction"]["evaluate"]["field"] =
    @benchmarkable eval_vec_batch($recon_field, $recon_points)
SUITE["reconstruction"]["evaluate"]["poloidal_flux"] =
    @benchmarkable eval_scalar_batch(poloidal_flux, $recon_field, $recon_points)
SUITE["reconstruction"]["evaluate"]["current_density"] =
    @benchmarkable eval_vec_batch(current_density, $recon_field, $recon_j_points)

# 3D field reconstruction of a tilted external coil from a 4 × 4 × 4 probe array
# and soft normal-field samples on the box surface.
SUITE["reconstruction3d"] = BenchmarkGroup()
recon3d_problem = let
    domain = (-0.5, 0.5, -0.5, 0.5, -0.5, 0.5)
    normal = SVector(0.3, 0.2, 1.0) / sqrt(1.13)
    coil = CurrentLoop(1.0, 1.0e4, SVector(0.2, -0.1, 2.0), normal)
    truth = CurrentLoopAnalytic(coil)
    to_matrix(points) = [p[i] for p in points, i in 1:3]
    lattice(r) = vec([SVector(x, y, z) for x in r, y in r, z in r])
    probes = lattice(range(-0.45, 0.45; length = 4))
    positions = to_matrix(probes)
    B = to_matrix(truth.(probes))
    surface_points, surface_normals = SVector{3, Float64}[], SVector{3, Float64}[]
    for axis in 1:3, side in (-0.5, 0.5), a in (-0.35, 0.0, 0.35), b in (-0.35, 0.0, 0.35)
        tangential = (a, b)
        push!(
            surface_points,
            SVector(ntuple(k -> k == axis ? side : tangential[k - (k > axis)], 3))
        )
        push!(surface_normals, SVector(ntuple(k -> k == axis ? sign(side) : 0.0, 3)))
    end
    boundary_values = map(zip(surface_points, surface_normals)) do (p, n)
        sum(n .* truth(p))
    end
    common = (;
        domain, uncertainties = fill(1.0e-5, size(B)),
        boundary_positions = to_matrix(surface_points),
        boundary_normals = to_matrix(surface_normals), boundary_values,
        boundary_uncertainties = fill(2.0e-5, length(boundary_values)),
        regularization = 1.0e-8,
    )
    centers = to_matrix(lattice(range(-0.5, 0.5; length = 4)))
    (; positions, B, common, centers)
end

SUITE["reconstruction3d"]["fit"] = BenchmarkGroup()
SUITE["reconstruction3d"]["fit"]["wendland_lsmr"] = @benchmarkable(
    reconstruct_field_3d(p.positions, p.B; length_scale = 0.5, p.common...),
    setup = (p = $recon3d_problem)
)
SUITE["reconstruction3d"]["fit"]["gaussian_lsmr"] = @benchmarkable(
    reconstruct_field_3d(
        p.positions, p.B; kernel = :gaussian, length_scale = 1.2,
        basis_centers = p.centers, regularization_grid = (7, 7, 7), p.common...
    ),
    setup = (p = $recon3d_problem)
)

recon3d_field = reconstruct_field_3d(
    recon3d_problem.positions, recon3d_problem.B; length_scale = 0.5,
    recon3d_problem.common...
)
recon3d_r = SVector(0.1, -0.2, 0.3)
recon3d_points = vec(
    [
        SVector(x, y, z) for x in (-0.35, 0.0, 0.35),
            y in (-0.35, 0.0, 0.35), z in (-0.35, 0.0, 0.35)
    ]
)

SUITE["reconstruction3d"]["evaluate"] = BenchmarkGroup()
SUITE["reconstruction3d"]["evaluate"]["field"] =
    @benchmarkable eval_vec_batch($recon3d_field, $recon3d_points)
SUITE["reconstruction3d"]["evaluate"]["current_density"] =
    @benchmarkable eval_vec_batch(current_density, $recon3d_field, $recon3d_points)

SUITE["reconstruction3d"]["assembly"] = BenchmarkGroup()
SUITE["reconstruction3d"]["assembly"]["operators"] =
    @benchmarkable Magnetostatics._operators_3d($recon3d_field, $recon3d_r)
recon3d_rows = reduce(
    vcat, [
        first(Magnetostatics._operators_3d(recon3d_field, SVector{3}(p)))
            for p in eachrow(recon3d_problem.positions)
    ]
)
recon3d_ncoefficients = length(recon3d_field.coefficients)
SUITE["reconstruction3d"]["assembly"]["matrix"] = @benchmarkable(
    Magnetostatics._matrix_3d($recon3d_rows, $recon3d_ncoefficients, Float64)
)

SUITE["grad_shafranov"] = BenchmarkGroup()
gsr_problem = let
    B0, L, epsilon, offset = 20.0e-9, 1.0e6, 0.3, 0.2e6
    island = FadeevIsland(B0, L, L, epsilon)
    times = collect(range(-40.0, 40.0; length = 121))
    xs = 50.0e3 .* times
    B = map(xs) do x
        field = island(SVector(x + π * L, 0.0, offset))
        SVector(field[1], field[3], 5.0e-9)
    end
    pressure = [
        0.1e-9 + B0^2 * (1 - epsilon^2) /
            (2Magnetostatics.μ₀ * (cosh(offset / L) - epsilon * cos(x / L))^2)
            for x in xs
    ]
    common = (;
        crossing_velocity = SVector(50.0e3, 0.0, 0.0),
        invariant_axis = SVector(0.0, 0.0, 1.0),
        y = collect(range(-0.08L, 0.08L; length = 61)),
        spatial_degree = 10, pressure_degree = 6, axial_degree = 0,
    )
    (; times, B = [b[k] for b in B, k in 1:3], pressure, common)
end
SUITE["grad_shafranov"]["fit"] = @benchmarkable(
    reconstruct_grad_shafranov(p.times, p.B, p.pressure; p.common...),
    setup = (p = $gsr_problem)
)
gsr_field = reconstruct_grad_shafranov(
    gsr_problem.times, gsr_problem.B, gsr_problem.pressure; gsr_problem.common...
)
gsr_points = vec(
    [
        SVector(x, y, 0.0) for x in range(-1.4e6, 1.4e6; length = 8),
            y in range(-0.07e6, 0.07e6; length = 8)
    ]
)
SUITE["grad_shafranov"]["field"] =
    @benchmarkable eval_vec_batch($gsr_field, $gsr_points)
SUITE["grad_shafranov"]["current_density"] =
    @benchmarkable eval_vec_batch(current_density, $gsr_field, $gsr_points)
SUITE["grad_shafranov"]["flux"] =
    @benchmarkable eval_scalar_batch(flux_function, $gsr_field, $gsr_points)
SUITE["grad_shafranov"]["pressure"] =
    @benchmarkable eval_scalar_batch(thermal_pressure, $gsr_field, $gsr_points)

SUITE["topology"] = BenchmarkGroup()
topo_ah_coil = AntiHelmholtzCoil(1.0, 10.0)
SUITE["topology"]["locate_null"] =
    @benchmarkable locate_magnetic_null($topo_ah_coil, SA[0.08, -0.05, 0.06])
