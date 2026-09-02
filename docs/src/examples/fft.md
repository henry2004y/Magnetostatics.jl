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

Visualizing the 3D magnetic field:

```@example fft
using CairoMakie, UniformStreamlines, FastInterpolations

fig = Figure(size = (800, 750), fontsize=20)
ax = Axis3(fig[1, 1];
    xlabel="x", ylabel="y", zlabel="z", aspect=:data, title="FFT Solver Result (3D Streamlines)")

# Build interpolants with FillExtrap so that NaN separator points in the
# streamline paths do not trigger a DomainError inside colorize.
Bx, By, Bz = @view(B[1, :, :, :]), @view(B[2, :, :, :]), @view(B[3, :, :, :])
Bx_itp = linear_interp((x, y, z), Bx; extrap = FillExtrap(0.0))
By_itp = linear_interp((x, y, z), By; extrap = FillExtrap(0.0))
Bz_itp = linear_interp((x, y, z), Bz; extrap = FillExtrap(0.0))

# Use evenstream for uniform 3D field line placement using grid data directly
str = evenstream(x, y, z, Bx_itp, By_itp, Bz_itp; min_density=0.5, max_density=1.0)
c = colorize(str, :norm) 
streamlines!(ax, str; color=c, colormap=:turbo, linewidth = 1.5, with_arrows=true, arrows_spacing=0.4)

fig
```
