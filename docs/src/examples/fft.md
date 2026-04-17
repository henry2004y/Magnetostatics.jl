# [FFT Solver](@id fft_example)

Solve for **B** given a current density **J** on a uniform grid.

```@example fft
using Magnetostatics, StaticArrays, LinearAlgebra, FFTW
using CairoMakie

# Define Grid
Nx, Ny, Nz = 8, 8, 8
dx = 0.1
J = zeros(Float64, 3, Nx, Ny, Nz)

# Populate J with a wire along z at the center
# Use a Gaussian distribution to minimize aliasing
width = 2.5 * dx
I = 1.0 # Current [A]

x = (0:Nx-1) .* dx .- (Nx*dx/2)
y = (0:Ny-1) .* dx .- (Ny*dx/2)
z = (0:Nz-1) .* dx .- (Nz*dx/2)

set_current_wire!(J, x, y, z, SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 1.0), I, width)

# Solve
solver = FFTSolver()
B = solve(solver, J, dx)
println("Max B magnitude: ", maximum(sqrt.(sum(B.^2, dims=1))))
```

Visualizing the 3D magnetic field:

```@example fft
using UniformStreamlines

fig = Figure(size = (800, 750), fontsize=20)
ax = Axis3(fig[1, 1];
    xlabel="x", ylabel="y", zlabel="z", aspect=:data, title="FFT Solver Result (3D Streamlines)")

# Use evenstream for uniform 3D field line placement using grid data directly
str = evenstream(x, y, z, @view(B[1,:,:,:]), @view(B[2,:,:,:]), @view(B[3,:,:,:]); min_density=0.5, max_density=1.0)
c = colorize(str, :norm) 
streamlines!(ax, str; color=c, colormap=:turbo, linewidth = 1.5, with_arrows=true, arrows_spacing=0.4)

fig
```
