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
