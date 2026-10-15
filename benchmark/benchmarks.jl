using BenchmarkTools
using LinearSolve: QRFactorization
using Magnetostatics
using StaticArrays

const SUITE = BenchmarkGroup()

SUITE["solvers"] = BenchmarkGroup()

# Biot-Savart
SUITE["solvers"]["biot_savart"] = BenchmarkGroup()
solver, wire, r = let
    r = SVector(0.5, 0.5, 0.0)
    loop = CurrentLoop(1.0, 1.0, SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 1.0))
    # creates a 100-segment wire
    wire = discretize_loop(loop, 100)
    solver = BiotSavart()

    solver, wire, r
end

SUITE["solvers"]["biot_savart"]["discretized_loop"] = @benchmarkable solve($solver, $wire, $r)

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

# Vector Potential
SUITE["solvers"]["potential"] = BenchmarkGroup()
pot_solver, pot_r, pot_loop = let
    solver = VectorPotential()
    r = SVector(0.5, 0.0, 0.1)
    loop = CurrentLoop(1.0, 1.0, SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 1.0))
    solver, r, loop
end

SUITE["solvers"]["potential"]["current_loop"] = @benchmarkable solve($pot_solver, $pot_loop, $pot_r)

# Analytical
SUITE["solvers"]["analytical"] = BenchmarkGroup()
r, loop = let
    r = SVector(0.5, 0.0, 0.1)
    # radius=3, current=1
    loop = CurrentLoop(3.0, 1.0, SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 1.0))

    r, loop
end

SUITE["solvers"]["analytical"]["current_loop"] = @benchmarkable getB_loop($r, $loop)

r_dipole, dipole = let
    r = SVector(0.5, 0.0, 0.1)
    dipole = Dipole(SVector(0.0, 0.0, 1.0))
    r, dipole
end

SUITE["solvers"]["analytical"]["dipole"] = @benchmarkable ($dipole)($r_dipole)

# Poisson Solver
SUITE["solvers"]["poisson"] = BenchmarkGroup()
poisson_solver, poisson_J = let
    xs = range(0.0, 1.0; length = 16)
    ys = range(0.0, 1.0; length = 16)
    zs = range(0.0, 1.0; length = 16)
    solver = PoissonSolver(xs, ys, zs)
    J = rand(3, 16, 16, 16)
    solver, J
end

SUITE["solvers"]["poisson"]["small_grid"] = @benchmarkable solve($poisson_solver, $poisson_J)

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
recon_r = SVector(0.9, 0.1, 0.05)

SUITE["reconstruction"]["evaluate"] = BenchmarkGroup()
SUITE["reconstruction"]["evaluate"]["field"] = @benchmarkable ($recon_field)($recon_r)
SUITE["reconstruction"]["evaluate"]["poloidal_flux"] =
    @benchmarkable poloidal_flux($recon_field, $recon_r)
SUITE["reconstruction"]["evaluate"]["current_density"] =
    @benchmarkable current_density($recon_field, $recon_r)

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

SUITE["reconstruction3d"]["evaluate"] = BenchmarkGroup()
SUITE["reconstruction3d"]["evaluate"]["field"] = @benchmarkable ($recon3d_field)($recon3d_r)
SUITE["reconstruction3d"]["evaluate"]["current_density"] =
    @benchmarkable current_density($recon3d_field, $recon3d_r)
