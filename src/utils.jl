# Utility functions

"""
    discretize_loop(radius, n_segments, current; center=[0,0,0], normal=[0,0,1])

Discretize a circular current loop into a `Wire` object.
"""
function discretize_loop(
        radius, n_segments, current;
        center = SVector(0.0, 0.0, 0.0), normal = SVector(0.0, 0.0, 1.0)
    )
    # Helper to create a loop
    θ = range(0, 2π, length = n_segments + 1)

    # Define rotation from z-axis to normal
    z = SVector(0.0, 0.0, 1.0)
    n = normalize(normal)

    # Calculate rotation matrix
    if n ≈ z
        R = I
    elseif n ≈ -z
        R = @SMatrix [1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 -1.0]
    else
        axis = cross(z, n)
        angle = acos(clamp(dot(z, n), -1.0, 1.0))
        u = normalize(axis)
        K = @SMatrix [
            0.0   -u[3]  u[2];
            u[3]   0.0  -u[1];
            -u[2]   u[1]  0.0
        ]
        R = I + sin(angle) * K + (1 - cos(angle)) * K * K
    end

    # Apply rotation and translation
    points = [R * SVector(radius * cos(t), radius * sin(t), 0.0) + center for t in θ]
    return Wire(points, current)
end

"""
    discretize_loop(loop::CurrentLoop, n_segments)

Discretize a `CurrentLoop` object into a `Wire` object.
"""
function discretize_loop(loop::CurrentLoop, n_segments)
    return discretize_loop(
        loop.radius, n_segments, loop.current;
        center = loop.center, normal = loop.normal
    )
end

"""
    discretize_knot(knot::TorusKnot, n_segments) -> Wire

Discretize a torus knot into a `Wire` object.
"""
function discretize_knot(knot::TorusKnot{T}, n_segments) where {T}
    # Create local coordinate system
    n = knot.normal
    u = knot.up
    v = cross(n, u)

    points = Vector{SVector{3, T}}(undef, n_segments + 1)
    for i in 0:n_segments
        t = 2π * i / n_segments
        # Parametric coordinates in local frame
        x_loc = (knot.R + knot.r * cos(knot.q * t)) * cos(knot.p * t)
        y_loc = (knot.R + knot.r * cos(knot.q * t)) * sin(knot.p * t)
        z_loc = knot.r * sin(knot.q * t)

        # Transform to global frame
        points[i + 1] = knot.center + x_loc * u + y_loc * v + z_loc * n
    end
    return Wire(points, knot.current)
end

"""
    Wire(knot::TorusKnot, n_segments=200)

Alternative constructor for `Wire` from a `TorusKnot`.
"""
function Wire(knot::TorusKnot, n_segments = 200)
    return discretize_knot(knot, n_segments)
end

"""
    set_current_wire!(J, x, y, z, point, direction, current, width)

Set the current density contribution of a straight wire with a Gaussian profile to `J`.

# Arguments
- `J`: 4D array of size (3, Nx, Ny, Nz) representing the current density components.
- `x`, `y`, `z`: 1D arrays identifying the grid coordinates.
- `point`: A point on the wire (Cartesian coordinates).
- `direction`: Unit vector of the wire direction.
- `current`: Total current in [A].
- `width`: Gaussian width parameter.
"""
function set_current_wire!(J, x, y, z, point, direction, current, width)
    Nx, Ny, Nz = length(x), length(y), length(z)
    @assert size(J, 1) == 3 && size(J, 2) == Nx && size(J, 3) == Ny && size(J, 4) == Nz

    u = normalize(direction)
    p = SVector{3}(point)
    J0 = current / (π * width^2)

    inv_width_sq = 1 / width^2

    # Precompute constant factors for direction
    Jx_dir, Jy_dir, Jz_dir = u[1], u[2], u[3]

    # Precompute grid mesh if possible or just loop
    # Simple loop structure
    for k in 1:Nz
        zz = z[k]
        for j in 1:Ny
            yy = y[j]
            for i in 1:Nx
                xx = x[i]
                r = SVector(xx, yy, zz)

                # Vector from point p to r
                d_vec = r - p

                # Projection along the wire direction
                proj = dot(d_vec, u)

                # Perpendicular vector to the wire
                perp = d_vec - proj * u

                # Squared distance to the wire
                dist_sq = dot(perp, perp)

                val = exp(-dist_sq * inv_width_sq)
                factor = J0 * val

                J[1, i, j, k] += factor * Jx_dir
                J[2, i, j, k] += factor * Jy_dir
                J[3, i, j, k] += factor * Jz_dir
            end
        end
    end

    return
end

"""
    set_current_wire(x, y, z, point, direction, current, width)

Create a new current density array `J` for a straight wire with a Gaussian profile.
See `set_current_wire!`.
"""
function set_current_wire(x, y, z, point, direction, current, width)
    Nx, Ny, Nz = length(x), length(y), length(z)
    J = zeros(eltype(x), 3, Nx, Ny, Nz)
    set_current_wire!(J, x, y, z, point, direction, current, width)
    return J
end

"""
    sph2cart(r, θ, ϕ)

Convert from spherical to Cartesian coordinates vector.
"""
function sph2cart(r, θ, ϕ)
    sinθ, cosθ = sincos(θ)
    sinϕ, cosϕ = sincos(ϕ)
    return SVector{3}(r * sinθ * cosϕ, r * sinθ * sinϕ, r * cosθ)
end

@inline @inbounds sph2cart(x) = sph2cart(x[1], x[2], x[3])

"""
    dipole_fieldline(ϕ, L=2.5, nP=100)

Creates `nP` points on one field line of the magnetic field from a dipole. In a centered
dipole magnetic field model, the path along a given L shell can be described as r = L*cos²λ,
where r is the radial distance (in planetary radii) to a point on the line,
λ is its co-latitude, and L is the L-shell of interest.
"""
function dipole_fieldline(ϕ, L = 2.5, nP = 100)
    xyz = [sph2cart(L * sin(θ)^2, θ, ϕ) for θ in range(0, stop = π, length = nP)]
    x = getindex.(xyz, 1)
    y = getindex.(xyz, 2)
    z = getindex.(xyz, 3)
    return (x, y, z)
end

"""
    image_source(source::Wire, bc::ConductingWall) -> Wire

Return the image `Wire` of `source` reflected across the planar conducting wall
defined by `bc`.  The image current is negated so that the tangential magnetic
field vanishes at the wall surface.
"""
function image_source(source::Wire{T}, bc::ConductingWall) where {T}
    ax = bc.axis
    wall = T(bc.position)
    reflected = map(source.points) do p
        delta = p[ax] - wall
        coords = ntuple(i -> i == ax ? wall - delta : p[i], 3)
        SVector{3, T}(coords)
    end
    return Wire{T}(reflected, -source.current)
end

"""
    compute_curl!(B, A, xs, ys, zs)

Compute `B = ∇ × A` in-place using second-order central differences on a
uniform Cartesian grid with potentially anisotropic spacing.

Interior points are updated; boundary layers are left unchanged (zero).

# Arguments
- `B`: output 4-D array of size `(3, Nx, Ny, Nz)`.
- `A`: input 4-D array of size `(3, Nx, Ny, Nz)`.
- `xs`, `ys`, `zs`: grid coordinate ranges (spacing extracted via `step`).
"""
function compute_curl!(B, A, xs, ys, zs)
    @assert size(B) == size(A) "B and A must have the same size!"
    @assert size(B, 1) == size(A, 1) == 3 "B and A must have 3 components!"
    dx, dy, dz = step(xs), step(ys), step(zs)
    inv2dx, inv2dy, inv2dz = inv(2dx), inv(2dy), inv(2dz)
    _, Nx, Ny, Nz = size(A)
    for k in 2:(Nz - 1), j in 2:(Ny - 1), i in 2:(Nx - 1)
        B[1, i, j, k] = (A[3, i, j + 1, k] - A[3, i, j - 1, k]) * inv2dy -
            (A[2, i, j, k + 1] - A[2, i, j, k - 1]) * inv2dz
        B[2, i, j, k] = (A[1, i, j, k + 1] - A[1, i, j, k - 1]) * inv2dz -
            (A[3, i + 1, j, k] - A[3, i - 1, j, k]) * inv2dx
        B[3, i, j, k] = (A[2, i + 1, j, k] - A[2, i - 1, j, k]) * inv2dx -
            (A[1, i, j + 1, k] - A[1, i, j - 1, k]) * inv2dy
    end
    return B
end

"""
    compute_curl(A, xs, ys, zs)

Allocate and return `B = ∇ × A`.  See [`compute_curl!`](@ref) for details.
"""
function compute_curl(
        A::AbstractArray{T, 4},
        xs::AbstractRange, ys::AbstractRange, zs::AbstractRange
    ) where {T}
    B = zeros(T, size(A))
    return compute_curl!(B, A, xs, ys, zs)
end
