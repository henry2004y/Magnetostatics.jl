# Examples

## Analytical Fields

Current configurations with known analytical solutions are useful for testing and benchmarking.

```julia
using Magnetostatics, StaticArrays, LinearAlgebra, FFTW
```

### Magnetic Dipole
Calculate the field of a magnetic dipole moment $\mathbf{M} = (0,0,1)$ at various points.

```julia
M = SVector(0.0, 0.0, 1.0)
dipole = Dipole(M)

r = SVector(1.0, 0.0, 0.0) # Point on x-axis
B = dipole(r)
println("B at $r: $B")
```

### Current Loop
Calculate the field of a circular current loop.

```julia
# Define loop parameters
radius = 1.0
current = 10.0
center = SVector(0.0, 0.0, 0.0)
normal = SVector(0.0, 0.0, 1.0)

# Create loop object
loop = CurrentLoop(radius, current, center, normal)

# Query field
r = SVector(0.0, 0.0, 0.5)
B = getB_loop(r, loop)
```

### Magnetic Mirror
Two co-axial current loops separated by a distance.

```julia
distance = 2.0
a = 1.0 # Radius
I = 100.0 # Current * windings

# Field at the center (0,0,0)
B = getB_mirror(0.0, 0.0, 0.0, distance, a, I)
```

### Tokamak Coils
Magnetic field from a Tokamak topology consisting of 16 toroidal field coils and a plasma current.

```julia
a = 1.0  # Coil radius
b = 2.0  # Major radius offset
ICoils = 1000.0
IPlasma = 500.0

# Query point
x, y, z = 3.0, 0.0, 0.0
B = getB_tokamak_coil(x, y, z, a, b, ICoils, IPlasma)
```

### Harris Sheet
A current sheet model often used in space physics ($B_x(z) = B_0 \tanh(z/L)$).

```julia
B0 = 1.0       # Asymptotic field strength
L = 2.0        # Half-width
sheet = HarrisSheet(B0, L)

r = SVector(0.0, 0.0, 1.0)
B = sheet(r)
```

### Z-Pinch
Magnetic field from a current running through a straight wire of radius $a$.

```julia
I = 100.0 # Current [A]
a = 0.1   # Radius [m]

x, y, z = 0.2, 0.0, 0.0
B = getB_zpinch(x, y, z, I, a)
```

### Tokamak with q-profile
Reconstruct the magnetic field distribution from a safety factor ($q$) profile.

```julia
# Define q-profile as a function of normalized radius r/a
q_profile(r_norm) = 1.1 + r_norm^2

a = 1.0   # Minor radius
R0 = 3.0  # Major radius
B0 = 2.0  # Toroidal field on axis

# Query point inside plasma
x, y, z = 3.5, 0.0, 0.0
B = getB_tokamak_profile(x, y, z, q_profile, a, R0, B0)
```

## Numerical Solvers

### Biot-Savart Solver
For arbitrary wire geometries, discretize the path and sum the contributions.

```julia
# Define a circular loop and discretize it
loop = CurrentLoop(1.0, 1.0, [0, 0, 0], [0, 0, 1])
wire = discretize_loop(loop, 100)

# Define the solver
solver = BiotSavart()

# Solve for B at a point
r = SVector(0.0, 0.0, 0.5)
B = solve(solver, wire, r)
```

### FFT Solver
Solve for **B** given a current density **J** on a uniform grid.

```julia
# Define Grid
Nx, Ny, Nz = 64, 64, 64
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
```
