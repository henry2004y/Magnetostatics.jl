# [FFT Solver](@id fft_example)

Solve for **B** given a current density **J** on a uniform grid.

```@example fft
using Magnetostatics, StaticArrays, FFTW

# Define Grid
Nx, Ny, Nz = 8, 8, 8
dx = 0.1
J = zeros(Float64, 3, Nx, Ny, Nz)

# Populate J with a wire along z at the center
# Use a Gaussian distribution to minimize aliasing
width = 2.5 * dx
I = 1.0 # Current [A]

x = LinRange(-Nx * dx / 2, Nx * dx / 2 - dx, Nx)
y = LinRange(-Ny * dx / 2, Ny * dx / 2 - dx, Ny)
z = LinRange(-Nz * dx / 2, Nz * dx / 2 - dx, Nz)

set_current_wire!(J, x, y, z, SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 1.0), I, width)

# Solve
solver = FFTSolver()
B = solve(solver, J, dx)
println("Max B magnitude: ", maximum(sqrt.(sum(B.^2, dims=1))))
```

Visualizing the 3D magnetic field with FieldTracer.jl's structured-grid API:
the component arrays are indexed as `(x, y, z)` (`gridtype = "ndgrid"`).
FieldTracer interpolates the grid data internally, so no separate interpolants
are needed. Unlike the callable-field API, this API returns three coordinate
vectors, not a result with `points` and `status`.
We select seeds explicitly on three planes and color each line by its seed radius.
The physical step length is `dx / 20`; a finite point budget traces approximately
one turn per seed, or stops earlier if the line leaves the grid.

```@example fft
using CairoMakie, FieldTracer

fig = Figure(; size = (800, 750), fontsize = 20)
ax = Axis3(fig[1, 1];
    xlabel = "x [m]", ylabel = "y [m]", zlabel = "z [m]", aspect = :data,
    title = "FFT Solver Result (FieldTracer Field Lines)")

Bx, By, Bz = @view(B[1, :, :, :]), @view(B[2, :, :, :]), @view(B[3, :, :, :])
ds = dx / 20
radii = (0.08, 0.16, 0.24)
for seed_z in (-0.2, 0.0, 0.2), radius in radii
    path = trace(
        Bx, By, Bz, radius, 0.0, seed_z, x, y, z;
        ds, maxstep = ceil(Int, 2π * radius / ds) + 1,
        gridtype = "ndgrid", direction = "forward"
    )
    @assert length(first(path)) > 1
    @assert all(all(isfinite, coordinates) for coordinates in path)
    lines!(
        ax, path...; color = fill(radius, length(first(path))),
        colorrange = extrema(radii), colormap = :turbo, linewidth = 1.5
    )
end
Colorbar(fig[1, 2]; limits = extrema(radii), colormap = :turbo, label = "Seed radius [m]")

fig
```
