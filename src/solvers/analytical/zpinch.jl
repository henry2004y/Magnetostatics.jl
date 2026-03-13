"""
    getB_zpinch(x, y, z, I, a) -> SVector{3}
    getB_zpinch(r, I, a) -> SVector{3}

Get magnetic field from a Z-pinch configuration.
Reference: [Z-pinch](https://en.wikipedia.org/wiki/Z-pinch)

# Arguments

  - `r`: location, vector of length 3 [m]
  - `x,y,z`: location [m]
  - `I::Float`: current in the wire [A].
  - `a::Float`: radius of the wire [m].
"""
@inline function getB_zpinch(x, y, z, I, a)
    return getB_zpinch(SVector(x, y, z), I, a)
end

@inline function getB_zpinch(pos, I, a)
    x, y = pos[1], pos[2]
    r = hypot(x, y)
    if r < a
        factor = μ₀ * I / (2π * a^2)
        Bx = -factor * y
        By = factor * x
    else
        factor = μ₀ * I / (2π * r^2)
        Bx = -factor * y
        By = factor * x
    end

    return SVector(Bx, By, 0.0)
end
