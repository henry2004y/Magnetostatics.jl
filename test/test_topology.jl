module test_topology

using Test
using Magnetostatics
using StaticArrays
using LinearAlgebra

struct LinearTestField <: AbstractMagneticField
    J::SMatrix{3, 3, Float64, 9}
    x0::SVector{3, Float64}
end

(f::LinearTestField)(r) = f.J * (SVector{3, Float64}(r) - f.x0)

@testset "Null Classification Taxonomy" begin
    @test classify_magnetic_null(@SMatrix [-1.0 0 0; 0 -1.0 0; 0 0 2.0]) === :A
    @test classify_magnetic_null(@SMatrix [1.0 0 0; 0 1.0 0; 0 0 -2.0]) === :B
    @test classify_magnetic_null(@SMatrix [-1.0 -2.0 0; 2.0 -1.0 0; 0 0 2.0]) === :As
    @test classify_magnetic_null(@SMatrix [1.0 -2.0 0; 2.0 1.0 0; 0 0 -2.0]) === :Bs
    @test classify_magnetic_null(@SMatrix [0.0 1.0 0; 1.0 0.0 0; 0 0 0.0]) === :X
    @test classify_magnetic_null(@SMatrix [0.0 -1.0 0; 1.0 0.0 0; 0 0 0.0]) === :O
    @test classify_magnetic_null(zeros(3, 3)) === :degenerate
    @test classify_magnetic_null(@SMatrix [1.0 0 0; 0 1.0 0; 0 0 1.0]) === :degenerate
end

@testset "3D Null Locator and Jacobian" begin
    let ah_b = AntiHelmholtzCoil(1.0, 10.0),
            ah_a = AntiHelmholtzCoil(1.0, -10.0),
            pt = SVector(0.12, -0.08, 0.15)
        jac_alloc(f, r) = (field_jacobian(f, r); @allocated field_jacobian(f, r))
        @test jac_alloc(ah_b, pt) == 0

        res_b = locate_magnetic_null(
            ah_b, SA[0.08, -0.05, 0.06];
            bounds = (-0.5, 0.5, -0.5, 0.5, -0.5, 0.5)
        )
        @test res_b.status === :converged
        @test res_b.kind === :B
        @test norm(res_b.point) < 1.0e-10
        @test res_b.residual < 1.0e-11

        res_a = locate_magnetic_null(ah_a, SA[-0.05, 0.07, -0.06])
        @test res_a.status === :converged
        @test res_a.kind === :A
        @test norm(res_a.point) < 1.0e-10

        # Transformed 3D null
        offset = SVector(0.3, -0.2, 0.15)
        R_y = @SMatrix [cos(0.4) 0.0 sin(0.4); 0.0 1.0 0.0; -sin(0.4) 0.0 cos(0.4)]
        ah_tr = translate_field(rotate_field(ah_a, R_y), offset)
        res_tr = locate_magnetic_null(ah_tr, offset + SA[0.04, 0.03, -0.05])
        @test res_tr.status === :converged
        @test res_tr.kind === :A
        @test res_tr.point ≈ offset atol = 1.0e-10

        # Spiral 3D null
        spiral = LinearTestField(
            @SMatrix([-1.0 -2.5 0.0; 2.5 -1.0 0.0; 0.0 0.0 2.0]),
            SA[0.2, -0.1, 0.3]
        )
        res_sp = locate_magnetic_null(spiral, SA[0.0, 0.0, 0.0])
        @test res_sp.status === :converged
        @test res_sp.kind === :As
        @test res_sp.point ≈ SA[0.2, -0.1, 0.3] atol = 1.0e-11

        # 2D null line from FadeevIsland (both X at x=0 and O at x=π*Lx)
        fad = FadeevIsland(1.0, 0.5, 1.0, 0.4)
        res_x = locate_magnetic_null(fad, SA[0.05, 0.1, 0.04])
        @test res_x.status === :converged
        @test res_x.kind === :X
        @test abs(res_x.point[1]) < 1.0e-10
        @test abs(res_x.point[3]) < 1.0e-10

        res_o = locate_magnetic_null(fad, SA[π - 0.08, 0.0, 0.04])
        @test res_o.status === :converged
        @test res_o.kind === :O
        @test abs(res_o.point[1] - π) < 1.0e-10
        @test abs(res_o.point[3]) < 1.0e-10

        @test locate_magnetic_null(NullField(), SA[0.0, 0.0, 0.0]).status === :degenerate
    end
end

@testset "Poloidal Critical Point Locator" begin
    let sol = SolovevEquilibrium(2.0, 1.5, 1.5, 1.2)
        res_o = locate_poloidal_critical_point(
            sol, (1.8, 0.1); bounds = (0.5, 3.0, -1.0, 1.0)
        )
        @test res_o.status === :converged
        @test res_o.kind === :O_point
        @test res_o.R ≈ 1.5 atol = 1.0e-10
        @test abs(res_o.Z) < 1.0e-10
    end

    let flare = FLARECoils(; I_pf = 45.0e3, I_ef = 12.0e3, I_p = 15.0e3)
        res_x = locate_poloidal_critical_point(
            flare, SVector(0.85, 0.0, 0.05); bounds = (0.3, 1.5, -0.5, 0.5)
        )
        @test res_x.status === :converged
        @test res_x.kind === :X_point
        @test res_x.residual < 1.0e-11
    end
end

end # module test_topology
