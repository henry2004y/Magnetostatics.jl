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
