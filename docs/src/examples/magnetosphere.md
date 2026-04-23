# Analytical Earth Magnetosphere Model

This example demonstrates how to build a complex analytical model of Earth's magnetosphere by superposing simpler magnetic field components. This approach is widely used in space physics to study the large-scale topology of the magnetosphere and phenomena like magnetic reconnection.

The total magnetic field is represented as the linear superposition of four distinct contributions:
$$\mathbf{B}_{\text{total}} = \mathbf{B}_{\text{int}} + \mathbf{B}_{\text{mp}} + \mathbf{B}_{\text{tail}} + \mathbf{B}_{\text{IMF}}$$

## 1. The Internal Field: The Planetary Dipole

Earth's intrinsic magnetic field is primarily dipolar. We use the `Dipole` model from `Magnetostatics.jl` with Earth's magnetic moment. In normalized units (Earth radii $R_E$), the field at the magnetic equator on the surface is approximately 31,000 nT.

## 2. The Magnetopause: Image Dipole Method

To model the shielding effect of the Chapman-Ferraro currents, we use a smooth scalar weighting function $w(\mathbf{r})$ that transitions from 1 inside the magnetosphere to 0 in the solar wind. Instead of applying this to the magnetic field directly, we apply it to the **magnetic vector potential** $\mathbf{A}$.

The total vector potential is:
$$\mathbf{A}_{\text{total}} = w(\mathbf{r}) \mathbf{A}_{\text{inner}} + (1 - w(\mathbf{r})) \mathbf{A}_{\text{outer}}$$
where $\mathbf{A}_{\text{inner}}$ contains the dipole, image dipole, and magnetotail, and $\mathbf{A}_{\text{outer}}$ represents the IMF. Taking the curl analytically ensures the total field is perfectly divergence-free ($\nabla \cdot \mathbf{B} = 0$):
$$\mathbf{B}_{\text{total}} = w \mathbf{B}_{\text{inner}} + (1 - w) \mathbf{B}_{\text{outer}} + \nabla w \times (\mathbf{A}_{\text{inner}} - \mathbf{A}_{\text{outer}})$$
The third term $\nabla w \times (\mathbf{A}_{\text{inner}} - \mathbf{A}_{\text{outer}})$ analytically generates the magnetopause surface currents within a transition layer of thickness $d$. We define the boundary as a paraboloid of revolution: $f(\mathbf{r}) = x - R_{\text{mp}} + (y^2 + z^2) / (2R_{\text{mp}})$, with $w(\mathbf{r}) = \frac{1}{2}(1 - \tanh(f/d))$.

## 3. The Magnetotail: Stretching and Reconnection

The Earth's magnetotail is an elongated region where magnetic field lines are stretched far into the nightside. We model this complex structure using two primary components:

1. **The Harris Current Sheet**: This provides the characteristic magnetic field reversal across the equatorial plane:
   $$B_x(z) = B_0 \tanh(z/L)$$
   This component represents the effect of the cross-tail current and is responsible for the elongated "tail-like" appearance of the field lines.
2. **The Inner Southward Field**: A realistic magnetosphere features a magnetic null point (X-line) where reconnection can occur. Since the Earth's intrinsic field is northward ($B_z > 0$) on the nightside equator, an X-line can only form if there is a competing southward field contribution. We introduce a background southward field $\mathbf{B}_{\text{tail}, z}$ to represent the combined effect of the tail current sheet and partial IMF penetration. The position of the X-line is determined by the balance between this component and the planetary dipole.

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
const DMP = 1.0      # Magnetopause thickness [RE]
const B_imf = UniformField(SVector(0.0, 0.0, -10.0)) # Southward IMF [nT]
const B_tail_z = UniformField(SVector(0.0, 0.0, -5.0)) # Southward field inside magnetosphere [nT]

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
struct SuperposedEarthField{D, T, U, S} <: AbstractMagneticField
    dipole_intrinsic::D
    dipole_mp::D
    tail::T
    imf::U
    tail_z::U
    image_pos::SVector{3, S}
    r_mp::S
    d_mp::S
end


function (f::SuperposedEarthField)(r)
    # 1. Geometry and weighting function
    # Paraboloid: x - Rmp + (y^2 + z^2)/(2Rmp) = 0
    x, y, z = r[1], r[2], r[3]
    dist = x - f.r_mp + (y^2 + z^2) / (2 * f.r_mp)
    w = 0.5 * (1 - tanh(dist / f.d_mp))
    
    # Gradient of weighting function
    dw_ddist = -0.5 * sech(dist / f.d_mp)^2 / f.d_mp
    grad_dist = SVector(1.0, y / f.r_mp, z / f.r_mp)
    grad_w = dw_ddist * grad_dist
    
    # 2. Inner components (Magnetosphere)
    # We include a southward B_tail_z representing tail current contributions
    B_inner = f.dipole_intrinsic(r) + f.dipole_mp(r - f.image_pos) + f.tail(r) + f.tail_z(r)
    A_inner = vector_potential(f.dipole_intrinsic, r) + 
              vector_potential(f.dipole_mp, r - f.image_pos) + 
              vector_potential(f.tail, r) +
              vector_potential(f.tail_z, r)
    
    # 3. Outer components (Solar Wind)
    B_outer = f.imf(r)
    A_outer = vector_potential(f.imf, r)
    
    # 4. Total B = w*B_inner + (1-w)*B_outer + grad_w x (A_inner - A_outer)
    return w * B_inner + (1 - w) * B_outer + cross(grad_w, A_inner - A_outer)
end

# Instantiate the model
mag_model = SuperposedEarthField(dipole_intrinsic, dipole_mp, tail, B_imf, B_tail_z, image_pos, RMP, DMP)

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

function get_B_xz(p)
    B = mag_model(SVector(p[1], 0.0, p[2]))
    return Point2f(B[1], B[3])
end

Bmag = [norm(mag_model(SVector(x, 0.0, z))) for x in xs, z in zs]

fig = Figure(size = (850, 700), fontsize=20)
ax = Axis(fig[1, 1], 
    xlabel=L"X [$R_E$]", ylabel=L"Z [$R_E$]", 
    title="Analytical Earth Magnetosphere Model",
    aspect=DataAspect())

# Heatmap of field strength
hm = heatmap!(ax, xs, zs, log10.(Bmag .+ 1.0), 
    colormap=:turbo, colorrange=(0, 5))
Colorbar(fig[1, 2], hm, label=L"$\log_{10}$(|B| [nT])")

# Streamlines
str = evenstream(xs, zs, get_B_xz; min_density=2.0)
streamlines!(ax, str; color = :black, linewidth = 1.0, with_arrows = true)

# Plot Earth
poly!(ax, Circle(Point2f(0, 0), RE), color=:gray)

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

