"""
    getB_mirror(x, y, z, distance, a, I1) :: SVector{3}
    getB_mirror(r, distance, a, I1) :: SVector{3}

Get magnetic field at `[x, y, z]` from a magnetic mirror generated from two coils.

# Arguments

  - `r`: location, vector of length 3 [m]
  - `x,y,z`: location [m]
  - `distance`: distance between solenoids in [m].
  - `a`: radius of each side coil in [m].
  - `I1`: current in the solenoid times number of windings in side coils.
"""
@inline function getB_mirror(x, y, z, distance, a, I1)
    return getB_mirror(SVector(x, y, z), distance, a, I1)
end

@inline function getB_mirror(
        r, distance, a, I1; center = SVector(0.0, 0.0, 0.0),
        normal = SVector(0.0, 0.0, 1.0)
    )
    c1 = center - 0.5 * distance * normal
    c2 = center + 0.5 * distance * normal
    cl1 = CurrentLoop(a, I1, c1, normal)
    cl2 = CurrentLoop(a, I1, c2, normal)
    B1 = getB_loop(r, cl1)
    B2 = getB_loop(r, cl2)

    return B1 + B2
end

"""
    getB_bottle(x, y, z, distance, a, b, I1, I2) :: SVector{3}
    getB_bottle(r, distance, a, b, I1, I2) :: SVector{3}

Get magnetic field from a magnetic bottle.
Reference: [wiki](https://en.wikipedia.org/wiki/Magnetic_mirror#Magnetic_bottles)

# Arguments

  - `r`: location, vector of length 3 [m]
  - `x,y,z`: location [m]
  - `distance::Float`: distance between solenoids in [m].
  - `a::Float`: radius of each side coil in [m].
  - `b::Float`: radius of central coil in [m].
  - `I1::Float`: current in the solenoid times number of windings in side coils in [A].
  - `I2::Float`: current in the central solenoid times number of windings in the
    central loop in [A].
"""
@inline function getB_bottle(
        x, y, z, distance, a, b, I1, I2;
        center = SVector(0.0, 0.0, 0.0),
        normal = SVector(0.0, 0.0, 1.0)
    )
    return getB_bottle(SVector(x, y, z), distance, a, b, I1, I2; center, normal)
end

@inline function getB_bottle(
        r, distance, a, b, I1, I2; center = SVector(0.0, 0.0, 0.0),
        normal = SVector(0.0, 0.0, 1.0)
    )
    B = getB_mirror(r, distance, a, I1; center, normal)

    # Central loop
    cl3 = CurrentLoop(b, I2, center, normal)
    B3 = getB_loop(r, cl3)

    return B + B3
end
