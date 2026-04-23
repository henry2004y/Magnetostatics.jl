# Analytical Earth Magnetosphere Model

This example demonstrates how to build a complex analytical model of Earth's magnetosphere by superposing simpler magnetic field components. This approach is widely used in space physics to study the large-scale topology of the magnetosphere and phenomena like magnetic reconnection.

The total magnetic field is represented as the linear superposition of four distinct contributions:
$$\mathbf{B}_{\text{total}} = \mathbf{B}_{\text{int}} + \mathbf{B}_{\text{mp}} + \mathbf{B}_{\text{tail}} + \mathbf{B}_{\text{IMF}}$$

## 1. The Internal Field: The Planetary Dipole

Earth's intrinsic magnetic field is primarily dipolar. We use the `Dipole` model from `Magnetostatics.jl` with Earth's magnetic moment. In normalized units (Earth radii $R_E$), the field at the magnetic equator on the surface is approximately 31,000 nT.

## 2. The Magnetopause: Image Dipole Method

To model the shielding effect of the Chapman-Ferraro currents at the magnetopause, we employ the method of images. By placing an "image dipole" upstream in the solar wind, we can cancel the normal component of the magnetic field at a specified standoff distance $R_{\text{mp}}$. For a subsolar standoff distance $R_{\text{mp}}$ along the x-axis, the image dipole is placed at $x = 2R_{\text{mp}}$.

## 3. The Magnetotail: The Harris Current Sheet

The nightside magnetotail is modeled using the Harris current sheet, which provides a magnetic field reversal across the equatorial plane:
$$B_x(z) = B_0 \tanh(z/L)$$

## 4. The Interplanetary Magnetic Field (IMF)

The IMF is modeled as a uniform background magnetic field. Adding a southward IMF ($B_z < 0$) creates magnetic nulls and enables the study of magnetic reconnection.

## Implementation

```@example earth_mag
using Magnetostatics, StaticArrays, LinearAlgebra
using CairoMakie, UniformStreamlines, GeometryBasics

# Constants (normalized to RE = 1 and nT)
const RE = 1.0
const BEQ = 31000.0  # Field at equator [nT]
const RMP = 10.0     # Magnetopause standoff distance [RE]
const B0_tail = 20.0 # Tail field strength [nT]
const L_tail = 2.0   # Tail sheet half-thickness [RE]
const B_imf = SVector(0.0, 0.0, -10.0) # Southward IMF [nT]

# Magnetic moments
# We scale the moment M such that (μ0_4π * M / RE^3) = BEQ
# In Magnetostatics.jl, μ0_4π is 1e-7
M_scaling = Magnetostatics.μ0_4π
M_val = BEQ / M_scaling

# Internal dipole (Southward moment for Northward field at surface)
M_int = SVector(0.0, 0.0, -M_val * RE^3)
image_pos = SVector(2 * RMP, 0.0, 0.0)
M_mp = M_int # Image dipole moment for cancellation

# Define components
dipole_intrinsic = Dipole(M_int)
dipole_mp = Dipole(M_mp)
tail = HarrisSheet(B0_tail, L_tail)

"""
    SuperposedEarthField <: AbstractMagneticField

A composite field model for Earth's magnetosphere.
"""
struct SuperposedEarthField{D, T, S} <: AbstractMagneticField
    dipole_intrinsic::D
    dipole_mp::D
    tail::T
    B_imf::SVector{3, S}
    image_pos::SVector{3, S}
end
function (f::SuperposedEarthField)(r)
    B_int = f.dipole_intrinsic(r)
    B_mp = f.dipole_mp(r - f.image_pos)
    B_tail = f.tail(r)
    return B_int + B_mp + B_tail + f.B_imf
end

# Instantiate the model
mag_model = SuperposedEarthField(dipole_intrinsic, dipole_mp, tail, B_imf, image_pos)

# Query at a point (e.g., 5 RE on the nightside)
r_test = SVector(-5.0, 0.0, 0.0)
B = mag_model(r_test)
println("B at $r_test: $B [nT]")
```

## Visualization

We visualize the magnetosphere in the X-Z plane (GSM coordinates), showing the compressed dayside and the elongated nightside tail.

```@example earth_mag
xs = range(-30, 15, length=100)
zs = range(-20, 20, length=100)

@inbounds function get_B_xz(p)
    B = mag_model(SVector(p[1], 0.0, p[2]))
    return Point2f(B[1], B[3])
end

Bmag = [norm(mag_model(SVector(x, 0.0, z))) for x in xs, z in zs]

fig = Figure(size = (900, 600), fontsize=20)
ax = Axis(fig[1, 1], 
    xlabel="X [RE]", ylabel="Z [RE]", 
    title="Analytical Earth Magnetosphere Model",
    aspect=DataAspect())

# Heatmap of field strength
hm = heatmap!(ax, xs, zs, log10.(Bmag .+ 1.0), 
    colormap=:viridis, colorrange=(0, 5))
Colorbar(fig[1, 2], hm, label="log10(|B| [nT])")

# Streamlines
str = evenstream(xs, zs, get_B_xz; min_density=2.0)
streamlines!(ax, str; color = :white, linewidth = 1.0, with_arrows = true)

# Plot Earth
poly!(ax, Circle(Point2f(0, 0), RE), color=:blue)

limits!(ax, -30, 15, -20, 20)
fig
```

## 3D Field Lines

We can also visualize the field lines in 3D to see the global topology of the magnetosphere.

```@example earth_mag
# 3D Field Lines Visualization
fig = Figure(size = (800, 700), fontsize=20)
ax = Axis3(fig[1, 1], title="3D Magnetosphere Field Lines", 
    aspect=:data, azimuth=1.55π, elevation=0.03π)

# Define the grid for 3D streamlines
xs_3d = range(-25, 10, length=21)
ys_3d = range(-15, 15, length=21)
zs_3d = range(-15, 15, length=21)

# Streamlines (sampling the field in 3D)
str = evenstream(xs_3d, ys_3d, zs_3d, p -> mag_model(p); min_density=0.8, max_density=2.0)
streamlines!(ax, str; color = :gray, linewidth = 1.0)

# Add Earth
mesh!(ax, Sphere(Point3f(0), RE), color=:blue)

# Add magnetopause standoff point for context
scatter!(ax, [Point3f(RMP, 0, 0)], color=:red, markersize=10, label="Standoff")

fig
```

