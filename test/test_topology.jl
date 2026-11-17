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

struct NaNField <: AbstractMagneticField end

(::NaNField)(_) = SVector(NaN, NaN, NaN)

struct CubicHarmonicField <: AbstractMagneticField end

function (::CubicHarmonicField)(r)
    x, y = r[1], r[2]
    return SVector(x^2 - y^2 - 0.16, -2.0 * x * y, 0.0)
end

@testset "Null Classification Taxonomy" begin
    @test classify_magnetic_null(@SMatrix [-1.0 0 0; 0 -1.0 0; 0 0 2.0]) === :A
    @test classify_magnetic_null(@SMatrix [1.0 0 0; 0 1.0 0; 0 0 -2.0]) === :B
    @test classify_magnetic_null(@SMatrix [-1.0 -2.0 0; 2.0 -1.0 0; 0 0 2.0]) === :As
    @test classify_magnetic_null(@SMatrix [1.0 -2.0 0; 2.0 1.0 0; 0 0 -2.0]) === :Bs
    @test classify_magnetic_null(@SMatrix [0.0 1.0 0; 1.0 0.0 0; 0 0 0.0]) === :X
    @test classify_magnetic_null(@SMatrix [0.0 -1.0 0; 1.0 0.0 0; 0 0 0.0]) === :O
    @test classify_magnetic_null(
        @SMatrix [0.05 0.05 0.0; -0.05 0.05 0.0; 0.0 0.0 0.0]
    ) === :degenerate
    @test classify_magnetic_null(
        @SMatrix [0.0 0.0 1.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
    ) === :degenerate
    @test classify_magnetic_null(
        @SMatrix [0.1 1.0 0.0; -0.01 0.1 0.0; 0.0 0.0 0.0]
    ) === :degenerate
    @test classify_magnetic_null(fill(NaN, 3, 3)) === :degenerate
    @test classify_magnetic_null(zeros(3, 3)) === :degenerate
    @test classify_magnetic_null(@SMatrix [1.0 0 0; 0 1.0 0; 0 0 1.0]) === :degenerate
    @test_throws ArgumentError classify_magnetic_null(zeros(2, 2))
    @test_throws ArgumentError classify_magnetic_null(zeros(3, 3); tol = -1.0)
    @test_throws ArgumentError classify_magnetic_null(zeros(3, 3); tol = Inf)
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
        @test res_b.eigenvalues ≈
            SVector{3, ComplexF64}(eigvals(Matrix(res_b.jacobian))) atol = 1.0e-12

        function null_alloc(f, r)
            locate_magnetic_null(f, r)
            return @allocated locate_magnetic_null(f, r)
        end
        cls_alloc(J) = (classify_magnetic_null(J); @allocated classify_magnetic_null(J))
        @test null_alloc(ah_b, SA[0.08, -0.05, 0.06]) == 0
        @test cls_alloc(res_b.jacobian) == 0

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
        @test res_sp.eigenvalues ≈
            SVector{3, ComplexF64}(eigvals(Matrix(res_sp.jacobian))) atol = 1.0e-12

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
        res_nan = locate_magnetic_null(NaNField(), SA[0.0, 0.0, 0.0])
        @test res_nan.status === :degenerate
        @test all(isnan, res_nan.eigenvalues)

        res_max = locate_magnetic_null(ah_b, SA[0.3, -0.25, 0.2]; maxiters = 1)
        @test res_max.status === :maxiters
        @test res_max.kind === :unconverged

        blocked = LinearTestField(SMatrix{3, 3, Float64}(I), SA[1.0, 0.0, 0.0])
        res_blk = locate_magnetic_null(
            blocked, SA[0.0, 0.0, 0.0]; bounds = (-1.0, 0.0, -1.0, 1.0, -1.0, 1.0)
        )
        @test res_blk.status === :maxiters
        @test res_blk.kind === :unconverged

        @test_throws ArgumentError locate_magnetic_null(ah_b, SA[0.0, 0.0])
        @test_throws ArgumentError locate_magnetic_null(ah_b, SA[0.0, 0.0, 0.0]; h = 0.0)
        @test_throws ArgumentError locate_magnetic_null(
            ah_b, SA[0.0, 0.0, 0.0]; tol = 0.0
        )
        @test_throws ArgumentError locate_magnetic_null(
            ah_b, SA[0.0, 0.0, 0.0]; maxiters = 0
        )
        @test_throws ArgumentError locate_magnetic_null(
            ah_b, SA[0.0, 0.0, 0.0]; bounds = (0.5, -0.5, 0.0, 1.0, 0.0, 1.0)
        )
        @test_throws ArgumentError locate_magnetic_null(
            ah_b, SA[2.0, 0.0, 0.0]; bounds = (-0.5, 0.5, -0.5, 0.5, -0.5, 0.5)
        )
        @test_throws ArgumentError locate_magnetic_null(ah_b, SA[NaN, 0.0, 0.0])
    end

    let box = (-0.5, 0.5, -0.5, 0.5, -0.5, 0.5),
            cub = CubicHarmonicField(),
            pts = reduce(
            vcat, [
                reshape([x, y, z], 1, 3)
                    for x in range(-0.4, 0.4; length = 4)
                    for y in range(-0.4, 0.4; length = 4)
                    for z in range(-0.4, 0.4; length = 4)
            ]
        ),
            B_obs = reduce(
            vcat, [reshape(collect(cub(SVector{3}(p))), 1, 3) for p in eachrow(pts)]
        )
        fit_cub = reconstruct_field_3d(
            pts, B_obs; domain = box, model = :curl_free, vacuum_degree = 3,
            regularization = 1.0e-14
        )
        res_back = locate_magnetic_null(fit_cub, SA[0.05, 0.0, 0.0])
        @test res_back.status === :converged
        @test res_back.point[1] ≈ 0.4 atol = 1.0e-9
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

        res_phi = locate_poloidal_critical_point(sol, (1.7, -0.08); phi = 0.4)
        @test res_phi.status === :converged
        @test res_phi.point ≈ SVector(1.5 * cos(0.4), 1.5 * sin(0.4), 0.0) atol = 1.0e-9

        res_3d_phi = locate_poloidal_critical_point(
            sol, SVector(1.7, 0.0, 0.05); phi = 0.3
        )
        @test res_3d_phi.status === :converged
        @test res_3d_phi.point ≈
            SVector(1.5 * cos(0.3), 1.5 * sin(0.3), 0.0) atol = 1.0e-9

        res_deg = locate_poloidal_critical_point(
            UniformField(SA[0.0, 0.0, 1.0]), (1.0, 0.0)
        )
        @test res_deg.status === :degenerate
        @test res_deg.kind === :degenerate

        res_max = locate_poloidal_critical_point(sol, (2.2, 0.4); maxiters = 1)
        @test res_max.status === :maxiters
        @test res_max.kind === :unconverged

        res_blk = locate_poloidal_critical_point(
            sol, (1.4, 0.0); bounds = (0.5, 1.4, -1.0, 1.0)
        )
        @test res_blk.status === :maxiters
        @test res_blk.kind === :unconverged

        @test_throws ArgumentError locate_poloidal_critical_point(sol, (1.0,))
        @test_throws ArgumentError locate_poloidal_critical_point(sol, (1.0, 0.0); h = -1.0)
        @test_throws ArgumentError locate_poloidal_critical_point(
            sol, (1.0, 0.0); tol = 0.0
        )
        @test_throws ArgumentError locate_poloidal_critical_point(
            sol, (1.0, 0.0); maxiters = 0
        )
        @test_throws ArgumentError locate_poloidal_critical_point(
            sol, (1.0, 0.0); bounds = (-0.5, 2.0, -1.0, 1.0)
        )
        @test_throws ArgumentError locate_poloidal_critical_point(sol, (1.0e-7, 0.0))
        @test_throws ArgumentError locate_poloidal_critical_point(
            sol, (0.2, 0.0); bounds = (0.5, 2.0, -1.0, 1.0)
        )
    end

    let flare = FLARECoils(; I_pf = 45.0e3, I_ef = 12.0e3, I_p = 15.0e3)
        res_x = locate_poloidal_critical_point(
            flare, SVector(0.85, 0.0, 0.05); bounds = (0.3, 1.5, -0.5, 0.5)
        )
        @test res_x.status === :converged
        @test res_x.kind === :X_point
        @test res_x.residual < 1.0e-11
    end

    let dom = (0.8, 1.7, -0.8, 0.8),
            sol = SolovevEquilibrium(2.0, 1.5, 1.5, 1.2),
            pts_rz = reduce(
            vcat, [
                reshape([R, Z], 1, 2)
                    for R in range(0.85, 1.65; length = 6)
                    for Z in range(-0.6, 0.6; length = 5)
            ]
        ),
            B_rz = reduce(
            vcat, [
                reshape(collect(sol(SVector(pts_rz[i, 1], 0.0, pts_rz[i, 2]))), 1, 3)
                    for i in axes(pts_rz, 1)
            ]
        )
        fit_ax = reconstruct_axisymmetric_field(
            pts_rz, B_rz; domain = dom, length_scale = 0.6, regularization = 1.0e-12
        )
        res_ax = locate_poloidal_critical_point(fit_ax, (0.88, 0.0))
        @test res_ax.status === :converged
        @test res_ax.R ≈ 1.5 atol = 0.03
    end
end

@testset "3D Field Covariance Completeness" begin
    let box = (-0.5, 0.5, -0.5, 0.5, -0.5, 0.5),
            probes = reduce(
            vcat, [
                reshape([x, y, z], 1, 3)
                    for x in range(-0.4, 0.4; length = 4)
                    for y in range(-0.4, 0.4; length = 4)
                    for z in range(-0.4, 0.4; length = 4)
            ]
        ),
            s_template = HarrisSheet(1.0, 0.4),
            dip_ext = Dipole(SA[2.0e3, -1.0e3, 3.0e3], SA[0.1, -0.2, 1.8]),
            truth = dip_ext + 0.5 * s_template,
            B_obs = reduce(
            vcat, [reshape(collect(truth(SVector{3}(p))), 1, 3) for p in eachrow(probes)]
        ),
            fit_all = reconstruct_field_3d(
            probes, B_obs;
            domain = box, model = :div_free, length_scale = 0.55,
            vacuum_degree = 2, dipole_positions = [0.1 -0.2 1.8],
            source_fields = (s_template,), uncertainties = fill(1.0e-4, size(probes)),
            regularization = 1.0e-6, store_covariance = true
        ),
            pt = SA[0.12, -0.09, 0.15]
        cov_B = field_covariance(fit_all, pt)
        @test issymmetric(cov_B)
        @test minimum(eigvals(Symmetric(Matrix(cov_B)))) >= -1.0e-14
        @test_throws DomainError field_covariance(fit_all, SA[1.0, 0.0, 0.0])

        fit_nocov = reconstruct_field_3d(
            probes, B_obs;
            domain = box, model = :div_free, length_scale = 0.55
        )
        @test_throws ArgumentError field_covariance(fit_nocov, pt)
    end
end

end # module test_topology
