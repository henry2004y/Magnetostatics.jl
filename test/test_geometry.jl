module test_geometry

using Test
using Magnetostatics
using StaticArrays
using LinearAlgebra

@testset "Magnetic Geometry & Drifts" begin
    # Test b_unit
    @test b_unit(SVector(0.0, 0.0, 0.0)) == SVector(0.0, 0.0, 0.0)
    b_test = b_unit(SVector(0.0, 3.0, 4.0))
    @test b_test ≈ SVector(0.0, 0.6, 0.8)
    @test norm(b_test) ≈ 1.0

    # Test on analytical ZPinch field
    # Outside wire: B_theta = mu0 I / (2pi r).
    # |B| = mu0 I / (2pi r).
    # d|B|/dr = -mu0 I / (2pi r^2) = -|B| / r.
    # Curvature of circular field lines: kappa = -r_hat / r = -(x, y, 0)/r^2.
    I_val = 10.0
    a_val = 0.2
    pinch = ZPinch(I_val, a_val)

    r_obs = SVector(1.0, 0.0, 0.0)
    @test field_magnitude(pinch, r_obs) ≈ norm(pinch(r_obs))

    ∇B = grad_b(pinch, r_obs)
    B_mag = norm(pinch(r_obs))
    @test ∇B[1] ≈ -B_mag / 1.0 rtol = 1.0e-3
    @test abs(∇B[2]) < 1.0e-8 && abs(∇B[3]) < 1.0e-8

    κ = magnetic_curvature(pinch, r_obs)
    # Curvature must point towards the wire center: -x_hat / r = SVector(-1.0, 0.0, 0.0)
    @test κ ≈ SVector(-1.0, 0.0, 0.0) rtol = 1.0e-3

    # Test drift velocities
    drifts = drift_velocities(pinch, r_obs; q = 1.0, m = 1.0, v_perp = 1.0, v_par = 1.0)
    # b is along y. ∇B is along -x. b × ∇B is along +z.
    # b × κ is also along +z.
    @test drifts.v_gradB[3] > 0.0
    @test drifts.v_curv[3] > 0.0
    @test drifts.v_drift ≈ drifts.v_gradB + drifts.v_curv

    # Test mirror ratio and loss cone angle
    mirror = MagneticMirror(2.0, 1.0, 100.0)
    r_throat = SVector(0.0, 0.0, 1.0) # at coil plane
    r_center = SVector(0.0, 0.0, 0.0) # at midpoint
    Rm = mirror_ratio(mirror, r_center, r_throat)
    @test Rm > 1.0
    θ_loss = loss_cone_angle(mirror, r_center, r_throat)
    @test 0.0 < θ_loss < π / 2
    @test sin(θ_loss) ≈ 1 / sqrt(Rm)

    # Test magnetic helicity
    A_vec = SVector(1.0, 2.0, 3.0)
    B_vec = SVector(4.0, 5.0, 6.0)
    @test magnetic_helicity(A_vec, B_vec) == dot(A_vec, B_vec)
end

end # module test_geometry
