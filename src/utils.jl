# Utility functions

"""
    discretize_loop(radius, n_segments, current; center=[0,0,0], normal=[0,0,1])

Discretize a circular current loop into a `Wire` object.
"""
function discretize_loop(radius, n_segments, current; center = SVector(0.0, 0.0, 0.0), normal = SVector(0.0, 0.0, 1.0))
    # Helper to create a loop
    θ = range(0, 2π, length = n_segments + 1)
    points = [SVector(radius * cos(t), radius * sin(t), 0.0) + center for t in θ]
    # Simple rotation logic could be added here if normal is not z-axis
    return Wire(points, current)
end
