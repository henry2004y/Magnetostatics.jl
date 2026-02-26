# [Z-Pinch](@id z_pinch_example)

Magnetic field from an infinite straight wire of radius $a$ carrying a uniform current $I$ in the $z$-direction. The magnetic field in cylindrical coordinates $(r, \theta, z)$ is given by:
```math
\mathbf{B}(r) = \begin{cases} \frac{\mu_0 I r}{2\pi a^2} \hat{\theta} & r \le a \\ \frac{\mu_0 I}{2\pi r} \hat{\theta} & r > a \end{cases}
```

```@example zpinch
using Magnetostatics, StaticArrays, LinearAlgebra
using CairoMakie

I = 100.0 # Current [A]
a = 0.1   # Radius [m]

x, y, z = 0.2, 0.0, 0.0
B = getB_zpinch(x, y, z, I, a)
println("B at $(SVector(x,y,z)): $B [T]")
```

Visualizing the field in the xy-plane (perpendicular to wire):

```@example zpinch
xs = range(-0.5, 0.5, length=51)
ys = range(-0.5, 0.5, length=51)

function field_xy_zpinch(x, y)
    B = getB_zpinch(x, y, 0.0, I, a)
    return Point2f(B[1], B[2])
end

Bmag = [norm(getB_zpinch(x, y, 0.0, I, a)) for x in xs, y in ys]

fig = Figure(size = (700, 600), fontsize=20)
ax = Axis(fig[1, 1], xlabel="x", ylabel="y", aspect=DataAspect(), title="Z-Pinch Field")

hm = heatmap!(ax, xs, ys, log10.(Bmag .+ 1e-9), colormap=:plasma)
Colorbar(fig[1, 2], hm, label="log10(|B|)")

streamplot!(ax, field_xy_zpinch, xs[1]..xs[end], ys[1]..ys[end];
    arrow_size = 8, linewidth = 1.5)

# Draw the wire circle
arc!(ax, Point2f(0,0), a, -π, π, color=:red, linestyle=:dash)

fig
```
