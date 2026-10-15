# [Current Sheets](@id current_sheet_examples)

Various current sheet models used in space physics and plasma physics.

## [Harris Sheet](@id harris_sheet_example)

A current sheet model often used in space physics, described by the exact field:
```math
\mathbf{B}(z) = B_0 \tanh\left(\frac{z}{L}\right) \hat{x}
```
where ``B_0`` is the asymptotic magnetic field strength and ``L`` is the half-width of the current sheet.

```@example currentsheet
using Magnetostatics, StaticArrays, LinearAlgebra
using CairoMakie, UniformStreamlines

B0 = 1.0       # Asymptotic field strength
L = 2.0        # Half-width
sheet = HarrisSheet(B0, L)

r = SVector(0.0, 0.0, 1.0)
B = sheet(r)
println("B at $r: $B [T]")
```

```@example currentsheet
xs = range(-5, 5, length=51)
zs = range(-2, 2, length=21)

field(x, z) = sheet(SA[x, 0.0, z])[SA[1, 3]]

fig = Figure(size = (800, 350), fontsize=20)
ax = Axis(fig[1, 1];
    xlabel="x", ylabel="z", aspect=DataAspect(), 
    title="Harris Sheet Magnetic Field")

# Generate uniform field lines
str = evenstream(xs, zs, field; min_density=1.5, max_density=3.0)
# Color by Bx component to show the reversal
color = colorize(str, :u)
cmax = maximum(abs, filter(!isnan, color))
crange = (-cmax, cmax)
colormap = Reverse(:RdBu)

sl = streamlines!(ax, str; 
    color, colormap, 
    linewidth=1.8,
    with_arrows=true)

Colorbar(fig[1, 2]; limits=crange, colormap, label="Bx [T]")

fig
```

## [Asymmetric Harris Sheet](@id asymmetric_harris_sheet_example)

An asymmetric version of the Harris sheet where the magnetic field strength on either side can be different.

```math
B_x(z) = \frac{B_1 + B_2}{2} + \frac{B_1 - B_2}{2} \tanh\left(\frac{z}{L}\right)
```

```@example currentsheet
B1, B2, L = 2.0, 1.0, 2.0
sheet_asym = AsymmetricHarrisSheet(B1, B2, L)

zs = range(-5, 5, length=100)
Bxs = [sheet_asym(SVector(0.0, 0.0, z))[1] for z in zs]

fig = Figure()
ax = Axis(fig[1, 1], xlabel="z", ylabel="Bx", title="Asymmetric Harris Sheet Profile")
lines!(ax, zs, Bxs)
fig
```

## [Force-Free Harris Sheet](@id force_free_harris_sheet_example)

A force-free current sheet configuration where the magnetic pressure is constant.

```math
\begin{aligned}
B_x(z) &= B_0 \tanh(z/L) \\
B_y(z) &= B_0 \operatorname{sech}(z/L)
\end{aligned}
```

```@example currentsheet
B0, L = 1.0, 2.0
sheet_ff = ForceFreeHarrisSheet(B0, L)

zs = range(-5, 5, length=100)
Bxs = [sheet_ff(SVector(0.0, 0.0, z))[1] for z in zs]
Bys = [sheet_ff(SVector(0.0, 0.0, z))[2] for z in zs]
Bzs =  @. sqrt(Bxs^2 + Bys^2)

fig = Figure()
ax = Axis(fig[1, 1], xlabel="z", ylabel="B", title="Force-Free Harris Sheet")
lines!(ax, zs, Bxs, label="Bx")
lines!(ax, zs, Bys, label="By")
lines!(ax, zs, Bzs, label="|B|", linestyle=:dash)
axislegend(ax)
fig
```

## [Guide-Field Harris Sheet](@id guide_field_sheet_example)

In magnetic reconnection studies (such as at the Earth's magnetopause, in tokamak tearing modes, and in laboratory experiments), the reconnecting current sheet often carries an out-of-plane guide field ``B_g`` and a normal field ``B_n``:
```math
B_x(z) = B_0 \tanh\left(\frac{z}{L}\right), \quad B_y(z) = B_g, \quad B_z(z) = B_n
```
The corresponding current density is directed along the ``y``-axis:
```math
J_y(z) = \frac{B_0}{\mu_0 L} \operatorname{sech}^2\left(\frac{z}{L}\right)
```
Unlike the force-free sheet above, the uniform guide field does not compensate the drop of ``|B_x|`` at the center, so the total field strength still has a minimum of ``\sqrt{B_g^2 + B_n^2}`` there.

```@example currentsheet
B0 = 20.0e-9   # Asymptotic reconnecting field [T] (20 nT)
L = 500.0e3    # Sheet half-thickness [m] (500 km)
Bg = 10.0e-9   # Guide field [T] (10 nT)
Bn = 2.0e-9    # Normal field [T] (2 nT)

sheet_guide = GuideFieldHarrisSheet(B0, L; Bg, Bn)

r0 = SVector(0.0, 0.0, 0.0)
println("B at sheet center: ", sheet_guide(r0), " [T]")
println("J at sheet center: ", current_density(sheet_guide, r0), " [A/m²]")

zs = range(-3 * L, 3 * L, length=100)
Bs = [sheet_guide(SVector(0.0, 0.0, z)) .* 1e9 for z in zs]

fig = Figure(size=(750, 500), fontsize=18)
ax = Axis(fig[1, 1], xlabel="z / L", ylabel="Field Components [nT]",
    title="Guide-Field Harris Current Sheet Profiles")

lines!(ax, zs ./ L, getindex.(Bs, 1), label="Bx (reconnecting)", linewidth=3,
    color=:blue)
lines!(ax, zs ./ L, getindex.(Bs, 2), label="By (guide field)", linewidth=3,
    color=:red, linestyle=:dash)
lines!(ax, zs ./ L, getindex.(Bs, 3), label="Bz (normal)", linewidth=2,
    color=:green, linestyle=:dot)
lines!(ax, zs ./ L, norm.(Bs), label="|B|", linewidth=3, color=:black)

axislegend(ax, position=:rb)
fig
```

For tracing charged particle orbits (such as Speiser motion and meandering orbits) in reconnecting current sheets with guide and normal magnetic fields, see [TestParticle.jl](https://github.com/henry2004y/TestParticle.jl).

## [Bifurcated Harris Sheet](@id bifurcated_harris_sheet_example)

A current sheet model with two peaks in the current density, often observed in the Earth's magnetotail:
```math
B_x(z) = \frac{B_0}{2} \left[ \tanh\left(\frac{z - d}{L}\right) + \tanh\left(\frac{z + d}{L}\right) \right]
```

The bifurcation is most clearly seen in the current density profile, which exhibits two distinct peaks at ``z = \pm d`` when ``d`` is sufficiently large compared to ``L``.

```@example currentsheet
B0, L, d = 1.0, 1.0, 2.0
sheet_bif = BifurcatedHarrisSheet(B0, L, d)

zs = range(-6, 6, length=200)
Bxs = [sheet_bif(SVector(0.0, 0.0, z))[1] for z in zs]
# Current density jy is obtained directly from the model
jys = [current_density(sheet_bif, SVector(0.0, 0.0, z))[2] for z in zs]

fig = Figure(size = (800, 350))
ax1 = Axis(fig[1, 1], xlabel="z", ylabel="Bx [T]", title="Magnetic Field Profile")
lines!(ax1, zs, Bxs)

ax2 = Axis(fig[1, 2], xlabel="z", ylabel="jy [A/m²]", title="Current Density Profile")
lines!(ax2, zs, jys)

fig
```

## [Fadeev Island](@id fadeev_island_example)

The Fadeev island model describes a chain of magnetic islands (magnetic reconnection).

```math
A_y(x, z) = -B_0 L \ln [ \cosh(z/L) + \epsilon \cos(x/L_x) ]
```

```@example currentsheet
B0, L, Lx, ε = 1.0, 2.0, 10.0, 0.5
fadeev = FadeevIsland(B0, L, Lx, ε)

xs = range(-20, 20, length=100)
zs = range(-5, 5, length=50)

field(x, z) = fadeev(SA[x, 0.0, z])[SA[1, 3]]

fig = Figure(size = (800, 280), fontsize=20)
ax = Axis(fig[1, 1];
    xlabel="x", ylabel="z", aspect=DataAspect(),
    limits = (-20, 20, -5, 5), 
    title="Fadeev Island Magnetic Field Lines")

# Generate uniform field lines
str = evenstream(xs, zs, field; min_density=1.5, max_density=3.0)
# Color by Bx component
color = colorize(str, :u)
cmax = maximum(abs, filter(!isnan, color))
crange = (-cmax, cmax)
colormap = Reverse(:RdBu)

sl = streamlines!(ax, str; 
    color, colormap, 
    linewidth=1.5,
    with_arrows=true)

Colorbar(fig[1, 2]; limits=crange, colormap, label="Bx [T]")

fig
```
