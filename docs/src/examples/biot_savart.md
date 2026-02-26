# [Biot-Savart Solver](@id biot_savart_example)

The magnetic field for a finite wire segment is calculated using the algebraic form of the Biot-Savart law:
```math
\mathbf{B} = \frac{\mu_0 I}{4\pi} \frac{d\mathbf{l} \times \mathbf{a}}{|d\mathbf{l} \times \mathbf{a}|^2} \left( \frac{d\mathbf{l} \cdot \mathbf{a}}{|\mathbf{a}|} - \frac{d\mathbf{l} \cdot \mathbf{b}}{|\mathbf{b}|} \right)
```
where $\mathbf{a} = \mathbf{r} - \mathbf{r}_{start}$ and $\mathbf{b} = \mathbf{r} - \mathbf{r}_{end}$ are vectors from the ends of the segment to the observation point $\mathbf{r}$, and $d\mathbf{l} = \mathbf{r}_{end} - \mathbf{r}_{start}$ is the segment vector.

For arbitrary wire geometries, discretize the path and sum the contributions.

```@example biotsavart
using Magnetostatics, StaticArrays, LinearAlgebra
using CairoMakie

# Define a circular loop and discretize it
loop = CurrentLoop(1.0, 1.0, [0, 0, 0], [0, 0, 1])
wire = discretize_loop(loop, 100)

# Define the solver
solver = BiotSavart()

# Solve for B at a point
r = SVector(0.0, 0.0, 0.5)
B = solve(solver, wire, r)
println("B at $r: $B [T]")
```

Visualizing the result:

```@example biotsavart
xs = range(-2, 2, length=51)
zs = range(-2, 2, length=51)

function field_xz_bs(x, z)
    B = solve(solver, wire, SVector(x, 0.0, z))
    return Point2f(B[1], B[3])
end

fig = Figure(size = (700, 600), fontsize=20)
ax = Axis(fig[1, 1];
    xlabel="x", ylabel="z", aspect=DataAspect(), title="Biot-Savart Loop Field")

Bmag = [norm(solve(solver, wire, SVector(x, 0.0, z))) for x in xs, z in zs]
hm = heatmap!(ax, xs, zs, log10.(Bmag .+ 1e-9), colormap=:plasma)
Colorbar(fig[1, 2], hm, label="log10(|B|)")

streamplot!(ax, field_xz_bs, -2..2, -2..2;
    arrow_size = 8, linewidth = 1.5)

fig
```
