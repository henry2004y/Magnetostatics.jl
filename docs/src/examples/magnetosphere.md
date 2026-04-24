# Analytical Earth Magnetosphere Model

This example demonstrates how to build a complex analytical model of Earth's magnetosphere by superposing simpler magnetic field components. This approach is widely used in space physics to study the large-scale topology of the magnetosphere and phenomena like magnetic reconnection.

The total magnetic field is represented as the linear superposition of four distinct contributions:
$$\mathbf{B}_{\text{total}} = \mathbf{B}_{\text{int}} + \mathbf{B}_{\text{mp}} + \mathbf{B}_{\text{tail}} + \mathbf{B}_{\text{IMF}}$$

## 1. The Internal Field: The Planetary Dipole

Earth's intrinsic magnetic field is primarily dipolar. We use the `Dipole` model from `Magnetostatics.jl` with Earth's magnetic moment. In normalized units (Earth radii $R_E$), the field at the magnetic equator on the surface is approximately 31,000 nT.

## 2. The Magnetopause: Image Dipole Method

To model the shielding effect of the Chapman-Ferraro currents, we use a smooth scalar weighting function $w(\mathbf{r})$ that transitions from 1 inside the magnetosphere to 0 in the solar wind. Instead of applying this to the magnetic field directly, we apply it to the **magnetic vector potential** $\mathbf{A}$.

The total vector potential is:
$$
\mathbf{A}_{\text{total}} = w(\mathbf{r}) \mathbf{A}_{\text{inner}} + (1 - w(\mathbf{r})) \mathbf{A}_{\text{outer}}
$$
where $\mathbf{A}_{\text{inner}}$ contains the dipole, image dipole, and magnetotail, and $\mathbf{A}_{\text{outer}}$ represents the IMF. Taking the curl analytically ensures the total field is perfectly divergence-free ($\nabla \cdot \mathbf{B} = 0$):
$$
\mathbf{B}_{\text{total}} = w \mathbf{B}_{\text{inner}} + (1 - w) \mathbf{B}_{\text{outer}} + \nabla w \times (\mathbf{A}_{\text{inner}} - \mathbf{A}_{\text{outer}})
$$
The third term $\nabla w \times (\mathbf{A}_{\text{inner}} - \mathbf{A}_{\text{outer}})$ analytically generates the magnetopause surface currents within a transition layer of thickness $d$. We define the boundary as a paraboloid of revolution: $f(\mathbf{r}) = x - R_{\text{mp}} + (y^2 + z^2) / (2R_{\text{mp}})$, with $w(\mathbf{r}) = \frac{1}{2}(1 - \tanh(f/d))$.

## 3. The Magnetotail: Stretching and Reconnection

The Earth's magnetotail is an elongated region where magnetic field lines are stretched far into the nightside. We model this complex structure using two primary components:

1. **The Harris Current Sheet**: This provides the characteristic magnetic field reversal across the equatorial plane:
   $$B_x(z) = B_0 \tanh(z/L)$$
   This component represents the effect of the cross-tail current and is responsible for the elongated "tail-like" appearance of the field lines.
2. **The Inner Southward Field**: A realistic magnetosphere features a magnetic null point (X-line) where reconnection can occur. Since the Earth's intrinsic field is northward ($B_z > 0$) on the nightside equator, an X-line can only form if there is a competing southward field contribution. We introduce a background southward field $\mathbf{B}_{\text{tail}, z}$ to represent the combined effect of the tail current sheet and partial IMF penetration. The position of the X-line is determined by the balance between this component and the planetary dipole.

## 4. The Bow Shock and Magnetosheath

To implement a bow shock without breaking the $\nabla\cdot\mathbf{B}=0$ constraint required by TestParticle.jl, we extend the smooth blending method to encompass three distinct spatial domains: the inner magnetosphere, the magnetosheath, and the pristine upstream solar wind.

### Defining the Nested Boundaries
We utilize two geometric surfaces. The magnetopause is a paraboloid $f(\mathbf{r})=0$. Upstream of the magnetopause lies the bow shock. Because the solar wind flow is deflected and expanded around the obstacle, empirical models typically represent the bow shock as a hyperboloid of revolution:
$$g(x, y, z)=x-R_{\text{bs}}+\sqrt{A(y^2+z^2)+B^2}-B$$
Here, $R_{\text{bs}}$ is the subsolar bow shock standoff distance (where $R_{\text{bs}} > R_{\text{mp}}$), and $A$ and $B$ control the flaring angle of the shock surface.
From these two surfaces, we generate two independent, smooth step functions using a hyperbolic tangent to dictate the transitions:
* $w_{\text{mp}}(\mathbf{r})$: Transitions from $1$ (inside the magnetosphere) to $0$ (in the magnetosheath).
* $w_{\text{bs}}(\mathbf{r})$: Transitions from $1$ (inside the magnetosheath) to $0$ (in the upstream solar wind).

### The Three-Region Vector Potential
To ensure strict divergence-free physics across multiple boundaries, we define vector potentials for all three regions:
* $\mathbf{A}_{\text{inner}}$: The planetary dipole, image dipole, and tail currents.
* $\mathbf{A}_{\text{sw}}$: The pristine interplanetary magnetic field (IMF).
* $\mathbf{A}_{\text{sheath}}$: The compressed and draped magnetic field within the magnetosheath.

We then nest the interpolations. First, construct the field completely external to the magnetopause:
$$\mathbf{A}_{\text{outer}}=w_{\text{bs}}(\mathbf{r})\mathbf{A}_{\text{sheath}}+(1-w_{\text{bs}}(\mathbf{r}))\mathbf{A}_{\text{sw}}$$
Next, blend this outer field with the inner magnetosphere field:
$$\mathbf{A}_{\text{total}}=w_{\text{mp}}(\mathbf{r})\mathbf{A}_{\text{inner}}+(1-w_{\text{mp}}(\mathbf{r}))\mathbf{A}_{\text{outer}}$$

### Calculating the Total Field and Shock Currents
When taking the curl to find $\mathbf{B}_{\text{total}}=\nabla\times\mathbf{A}_{\text{total}}$, the resulting field smoothly transitions through all three regions. The nested product rule automatically generates two distinct boundary current layers:
$$\mathbf{B}_{\text{total}}=w_{\text{mp}}\mathbf{B}_{\text{inner}}+(1-w_{\text{mp}})[w_{\text{bs}}\mathbf{B}_{\text{sheath}}+(1-w_{\text{bs}})\mathbf{B}_{\text{sw}}]+\mathbf{B}_{\text{currents}}$$
The $\mathbf{B}_{\text{currents}}$ term isolates the analytical currents generated at the boundaries:
$$\mathbf{B}_{\text{currents}}=\nabla w_{\text{mp}}\times(\mathbf{A}_{\text{inner}}-\mathbf{A}_{\text{outer}})+(1-w_{\text{mp}})[\nabla w_{\text{bs}}\times(\mathbf{A}_{\text{sheath}}-\mathbf{A}_{\text{sw}})]$$
The first cross product represents the Chapman-Ferraro magnetopause currents. The second cross product represents the shock currents flowing on the bow shock surface. By tuning the steepness of the $w_{\text{bs}}$ weighting function, you control the physical thickness of the shock ramp.

### Formulating the Magnetosheath Field
To define an accurate $\mathbf{A}_{\text{sheath}}$, the field cannot simply be uniform, as the plasma slows down and the magnetic field lines drape around the magnetopause obstacle. To satisfy the macroscopic jump conditions across the bow shock, the solar wind magnetic field is compressed by a ratio $r_c$.

A practical analytical approach is to map the upstream uniform vector potential $\mathbf{A}_{\text{sw}}$ through a continuous coordinate transformation $\mathbf{r}' = \mathcal{T}(\mathbf{r})$ that mimics potential flow around the paraboloid obstacle:
$$\mathbf{A}_{\text{sheath}}(\mathbf{r})\approx r_c \mathbf{A}_{\text{sw}}(\mathcal{T}(\mathbf{r}))$$
For this implementation, we map the uniform vector potential using $\mathcal{T}(\mathbf{r}) = (f(\mathbf{r}), y, z)$, where $f(\mathbf{r})$ is the magnetopause paraboloid geometry. This analytically stretches and compresses the vector potential, natively capturing the magnetic draping effect and the pile-up of flux on the dayside magnetopause, all while strictly preserving $\nabla\cdot\mathbf{B}=0$ everywhere in the simulation domain.

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

# Bow shock parameters
const RBS = 13.0     # Bow shock subsolar standoff distance [RE]
const ABS = 0.04     # Bow shock flaring parameter
const BBS = 2.0      # Bow shock hyperboloid parameter
const DBS = 0.5      # Bow shock thickness [RE]
const RC = 3.0       # Magnetosheath compression ratio

const B_imf = UniformField(SVector(0.0, 0.0, -10.0)) # Southward IMF [nT]
const B_tail_z = UniformField(SVector(0.0, 0.0, -5.0)) # Southward field inside magnetosphere [nT]

# Magnetic moments
# We scale the moment M such that (μ0_4π * M / RE^3) = BEQ, where μ0_4π is 1e-7.
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
    dipole_intrinsic::D     # Internal planetary dipole
    dipole_mp::D            # Image dipole for magnetopause shielding
    tail::T                 # Harris current sheet for the magnetotail
    imf::U                  # Interplanetary magnetic field (IMF)
    tail_z::U               # Southward field contribution inside magnetosphere
    image_pos::SVector{3, S} # Position of the image dipole
    r_mp::S                 # Magnetopause standoff distance
    d_mp::S                 # Magnetopause transition layer thickness
    r_bs::S                 # Bow shock standoff distance
    a_bs::S                 # Bow shock flaring parameter
    b_bs::S                 # Bow shock hyperboloid parameter
    d_bs::S                 # Bow shock transition layer thickness
    r_c::S                  # Magnetosheath compression ratio
end


# Helper functions for the magnetosphere components
@inline inner_field(f::SuperposedEarthField{D, T, U, S}, r) where {D, T, U, S} = 
    f.dipole_intrinsic(r) + f.dipole_mp(r .- f.image_pos) + f.tail(r) + f.tail_z(r)

@inline inner_potential(f::SuperposedEarthField{D, T, U, S}, r) where {D, T, U, S} = 
    vector_potential(f.dipole_intrinsic, r) + vector_potential(f.dipole_mp, r .- f.image_pos) + 
    vector_potential(f.tail, r) + vector_potential(f.tail_z, r)

@inline function sheath_potential(f::SuperposedEarthField, r)
    dist_mp = r[1] - f.r_mp + (r[2]^2 + r[3]^2) / (2 * f.r_mp)
    r_prime = SVector(dist_mp, r[2], r[3])
    return f.r_c * vector_potential(f.imf, r_prime)
end

@inline function sheath_field(f::SuperposedEarthField, r)
    B_imf = f.imf.B0
    bx = B_imf[1] - 0.5 * (B_imf[2] * r[2] + B_imf[3] * r[3]) / f.r_mp
    return f.r_c * SVector(bx, B_imf[2], B_imf[3])
end

@inline function outer_potential(f::SuperposedEarthField, r)
    rad_term = sqrt(f.a_bs * (r[2]^2 + r[3]^2) + f.b_bs^2)
    dist_bs = r[1] - f.r_bs + rad_term - f.b_bs
    w_bs = 0.5 * (1 - tanh(dist_bs / f.d_bs))
    return w_bs * sheath_potential(f, r) + (1 - w_bs) * vector_potential(f.imf, r)
end

@inline function outer_field(f::SuperposedEarthField, r)
    rad_term = sqrt(f.a_bs * (r[2]^2 + r[3]^2) + f.b_bs^2)
    dist_bs = r[1] - f.r_bs + rad_term - f.b_bs
    w_bs = 0.5 * (1 - tanh(dist_bs / f.d_bs))
    
    dw_ddist = -0.5 * sech(dist_bs / f.d_bs)^2 / f.d_bs
    grad_dist = SVector(1.0, f.a_bs * r[2] / rad_term, f.a_bs * r[3] / rad_term)
    grad_w_bs = dw_ddist * grad_dist
    
    A_sh, A_sw = sheath_potential(f, r), vector_potential(f.imf, r)
    B_sh, B_sw = sheath_field(f, r), f.imf(r)
    
    return w_bs * B_sh + (1 - w_bs) * B_sw + cross(grad_w_bs, A_sh - A_sw)
end

@inline function Magnetostatics.vector_potential(f::SuperposedEarthField{D, T, U, S}, r) where {D, T, U, S}
    dist_mp = r[1] - f.r_mp + (r[2]^2 + r[3]^2) / (2 * f.r_mp)
    w_mp = 0.5 * (1 - tanh(dist_mp / f.d_mp))
    return w_mp * inner_potential(f, r) + (1 - w_mp) * outer_potential(f, r)
end

@inline function (f::SuperposedEarthField{D, T, U, S})(r) where {D, T, U, S}
    # Geometry and weighting function
    x, y, z = r[1], r[2], r[3]
    dist_mp = x - f.r_mp + (y^2 + z^2) / (2 * f.r_mp)
    w_mp = 0.5 * (1 - tanh(dist_mp / f.d_mp))
    
    # Gradient of weighting function
    dw_ddist = -0.5 * sech(dist_mp / f.d_mp)^2 / f.d_mp
    grad_dist = SVector(one(S), y / f.r_mp, z / f.r_mp)
    grad_w_mp = dw_ddist * grad_dist
    
    # Field and Potential evaluation
    B_in, B_out = inner_field(f, r), outer_field(f, r)
    A_in, A_out = inner_potential(f, r), outer_potential(f, r)
    
    # Total B
    return w_mp * B_in + (1 - w_mp) * B_out + cross(grad_w_mp, A_in - A_out)
end

# Instantiate the model
mag_model = SuperposedEarthField(dipole_intrinsic, dipole_mp, tail, B_imf, B_tail_z, image_pos, RMP, DMP, RBS, ABS, BBS, DBS, RC)

# Query at a point (e.g., 5 RE on the nightside)
r_test = SVector(-5.0, 0.0, 0.0)
B = mag_model(r_test)
println("B at $r_test: $B [nT]")
```

## Visualization

We visualize the magnetosphere in the X-Z plane (GSM coordinates), showing the compressed dayside and the elongated nightside tail.

```@example earth_mag
xs = range(-30, 20, length=100)
zs = range(-20, 20, length=100)

function get_B_xz(p)
    B = mag_model(SVector(p[1], 0.0, p[2]))
    return Point2f(B[1], B[3])
end

Bmag = [norm(mag_model(SVector(x, 0.0, z))) for x in xs, z in zs]

fig = Figure(size = (900, 700), fontsize=20)
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

limits!(ax, -30, 20, -20, 20)
fig
```

## 3D Field Lines

We can also visualize the field lines in 3D to see the global topology of the magnetosphere.

```@example earth_mag
fig = Figure(size = (800, 700), fontsize=20)
ax = Axis3(fig[1, 1], title="3D Magnetosphere Field Lines", 
    aspect=:data, azimuth=1.55π, elevation=0.03π)

# Define the grid for 3D streamlines
xs_3d = range(-20, 18, length=21)
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

