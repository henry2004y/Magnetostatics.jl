module test_forwarddiff

using Test
using Magnetostatics
using ForwardDiff
using StaticArrays
using LinearAlgebra

function _ad_curl(f, r)
    J = ForwardDiff.jacobian(f, r)
    return SVector(J[3, 2] - J[2, 3], J[1, 3] - J[3, 1], J[2, 1] - J[1, 2])
end

@testset "Issue #79 Dipole MWE" begin
    let dip = Dipole(SA[0.0, 0.0, 8.0e22]),
            r0 = SA[1.5e7, 0.0, 0.0]
        J = ForwardDiff.jacobian(dip, r0)
        @test size(J) == (3, 3)
        @test all(isfinite, J)
        @test abs(tr(J)) <= 1.0e-12 * opnorm(J, Inf)
        @test isapprox(J, field_jacobian(dip, r0; h = 1.0); rtol = 1.0e-8)

        JA = ForwardDiff.jacobian(r -> vector_potential(dip, r), r0)
        @test all(isfinite, JA)
        @test isapprox(_ad_curl(r -> vector_potential(dip, r), r0), dip(r0); rtol = 1.0e-12)
    end
end

@testset "Analytical Fields ForwardDiff" begin
    let r_gen = SVector(1.2, -0.5, 0.4),
            r_cyl = SVector(3.4, 0.6, 0.35),
            r_wire_in = SVector(0.15, -0.1, 0.2),
            r_wire_out = SVector(0.8, 0.3, -0.2),
            r_parker = SVector(2.0e11, 1.0e11, 0.8e11),
            r_mag = SVector(5.0, 1.0, -1.0),
            r_tok = SVector(1.1, 0.15, 0.1),
            R_rot = @SMatrix(
            [
                cos(0.4) 0.0 sin(0.4);
                0.0 1.0 0.0;
                -sin(0.4) 0.0 cos(0.4)
            ]
        )
        loop = CurrentLoop(1.0, 2.5, SA[0.0, 0.0, 0.1], SA[0.1, 0.2, 1.0])
        mag = AnalyticalMagnetosphere(;
            intrinsic_dipole = Dipole(SA[0.0, 0.0, 1.0e4]),
            imf_field = UniformField(SA[0.0, 0.0, -5.0e-9]),
            has_shock = true
        )
        tok_prof = TokamakProfile(nr -> nr^2 + 2 * nr + 0.5, 0.5, 1.0, 1.0)
        sol_thin = Solenoid(0.5, 1.2, 500.0)
        sol_thick = Solenoid(0.5, 1.0, 400.0; thickness = 0.2)
        sph = SphereMagnet(SA[1.0e5, -2.0e5, 3.0e5], 0.4, SA[0.2, -0.1, 0.3])
        cyl = CylinderMagnet(
            0.3, 0.8, 5.0e5; center = SA[0.1, 0.2, -0.1], normal = SA[0.0, 1.0, 1.0]
        )
        cub_rot = CuboidMagnet(
            SA[0.6, 1.0, 1.4], SA[2.0e5, -1.0e5, 4.0e5];
            center = SA[0.1, -0.2, 0.3], rotation = R_rot
        )
        cub_axis = CuboidMagnet(SA[0.6, 1.0, 1.4], SA[2.0e5, -1.0e5, 4.0e5])
        flare_g = FLARECoils(; I_ef = 20.0e3, I_p = 5.0e3, B_guide = 0.05)
        iw = InfiniteWire(15.0, 0.4, SA[0.0, 0.0, 0.0], SA[0.2, 0.1, 1.0])
        sp_off = ScrewPinch(25.0, 0.4, 0.8; center = SA[0.05, -0.05])

        fields = [
            (NullField(), r_gen, true, true, true),
            (UniformField(SA[0.2, -0.4, 0.7]), r_gen, true, true, true),
            (HarrisSheet(1.5, 0.6), r_gen, true, true, true),
            (AsymmetricHarrisSheet(2.0, -1.0, 0.8), r_gen, true, true, true),
            (ForceFreeHarrisSheet(1.5, 0.7), r_gen, true, true, true),
            (BifurcatedHarrisSheet(1.8, 0.6, 0.9), r_gen, true, true, true),
            (FadeevIsland(1.2, 0.8, 1.1, 0.35), r_gen, true, false, true),
            (GuideFieldHarrisSheet(1.4, 0.6; Bg = 0.5, Bn = 0.2), r_gen, true, true, true),
            (Dipole(SA[1.0, -2.0, 3.0], SA[0.1, 0.2, -0.1]), r_gen, true, true, true),
            (iw, r_wire_in, true, true, true),
            (iw, r_wire_out, true, true, true),
            (CurrentLoopAnalytic(loop), r_gen, true, true, true),
            (MagneticMirror(1.2, 0.8, 1.5), r_gen, true, true, true),
            (MagneticBottle(1.6, 0.7, 1.0, 2.0, 0.8), r_gen, true, true, true),
            (HelmholtzCoil(0.8, 2.0), r_gen, true, true, true),
            (AntiHelmholtzCoil(0.8, 2.0), r_gen, true, true, true),
            (sol_thin, r_gen, true, true, true),
            (sol_thick, SA[0.52, 0.0, 0.1], true, true, true),
            (sol_thick, SA[0.15, 0.1, 0.1], true, true, true),
            (sph, SA[0.25, -0.12, 0.34], true, true, true),
            (sph, SA[0.9, 0.3, -0.2], true, true, true),
            (cyl, SA[0.25, 0.3, 0.1], true, true, true),
            (cub_rot, SA[0.8, -0.6, 1.2], false, true, true),
            (cub_rot, SA[0.15, -0.22, 0.35], false, true, true),
            (cub_axis, SA[0.3, 0.2, 1.2], false, true, true),
            (ZPinch(20.0, 0.4; center = SA[0.1, -0.1]), r_wire_in, true, true, true),
            (ZPinch(20.0, 0.4; center = SA[0.1, -0.1]), r_wire_out, true, true, true),
            (sp_off, r_wire_in, true, true, true),
            (sp_off, r_wire_out, true, true, true),
            (TaylorState(1.2, 1.7; center = SA[0.05, 0.1]), r_gen, true, true, true),
            (SolovevEquilibrium(2.0, 3.0, 1.5, 1.2), r_cyl, true, true, true),
            (ParkerSpiral(5.0e-9, 1.496e11), r_parker, true, true, true),
            (GoldHoyleFluxTube(1.3, 0.9; center = SA[0.1, -0.2]), r_gen, true, true, true),
            (FLARECoils(; I_ef = 20.0e3, I_p = 5.0e3), r_gen, true, true, true),
            (flare_g, r_gen, true, true, true),
            (TokamakCoils(0.5, 1.0, 1.0e3, 0.0), r_gen, false, false, true),
            (TokamakCoils(0.5, 1.0, 1.0e3, 500.0), r_gen, false, false, false),
            (tok_prof, r_tok, false, false, true),
            (mag, r_mag, true, false, true),
        ]

        for (field, pt, has_vp, check_j, div_free) in fields
            JB = ForwardDiff.jacobian(field, pt)
            @test all(isfinite, JB)
            if div_free
                @test abs(tr(JB)) <= 1.0e-11 * max(1.0, opnorm(JB, Inf))
            end

            if has_vp
                JA = ForwardDiff.jacobian(r -> vector_potential(field, r), pt)
                @test all(isfinite, JA)
                @test isapprox(
                    _ad_curl(r -> vector_potential(field, r), pt),
                    field(pt);
                    rtol = 1.0e-10,
                    atol = 1.0e-12
                )
            end

            if check_j
                JJ = ForwardDiff.jacobian(r -> current_density(field, r), pt)
                @test all(isfinite, JJ)
                @test isapprox(
                    _ad_curl(field, pt),
                    Magnetostatics.μ₀ * current_density(field, pt);
                    rtol = 1.0e-10,
                    atol = 1.0e-11 * max(1.0, norm(field(pt)))
                )
            end
        end

        for f_flux in (
                SolovevEquilibrium(2.0, 3.0, 1.5, 1.2),
                sol_thin,
                sol_thick,
                cyl,
                flare_g,
            )
            pt = f_flux isa SolovevEquilibrium ? r_cyl : r_gen
            g = ForwardDiff.gradient(r -> poloidal_flux(f_flux, r), pt)
            @test all(isfinite, g)
        end
    end
end

@testset "Axis, Midplane, and Special Points" begin
    let r_z = SVector(0.0, 0.0, 0.25),
            r_mid = SVector(0.4, -0.3, 0.0),
            loop_ax = CurrentLoopAnalytic(
            CurrentLoop(1.0, 2.5, SA[0.0, 0.0, 0.0], SA[0.0, 0.0, 1.0])
        ),
            sol = Solenoid(0.5, 1.2, 400.0),
            sol_thick = Solenoid(0.5, 1.0, 400.0; thickness = 0.2),
            mir = MagneticMirror(1.2, 0.8, 1.5),
            bot = MagneticBottle(1.6, 0.7, 1.0, 2.0, 0.8),
            ah = AntiHelmholtzCoil(0.8, 2.0),
            ts = TaylorState(1.2, 1.7),
            gh = GoldHoyleFluxTube(1.3, 0.9),
            zp = ZPinch(20.0, 0.4),
            sp = ScrewPinch(25.0, 0.4, 0.8),
            iw = InfiniteWire(15.0, 0.4),
            ps = ParkerSpiral(5.0e-9, 1.496e11),
            sph = SphereMagnet(SA[1.0e5, -2.0e5, 3.0e5], 0.4),
            tok_c = TokamakCoils(0.5, 1.0, 1.0e3, 500.0),
            tok_p = TokamakProfile(nr -> nr^2 + 2 * nr + 0.5, 0.5, 1.0, 1.0),
            fad = FadeevIsland(1.2, 0.8, 1.1, 0.35)

        for (field, pt) in (
                (loop_ax, r_z),
                (loop_ax, r_mid),
                (sol, r_z),
                (sol, r_mid),
                (sol_thick, r_z),
                (mir, r_z),
                (bot, r_z),
                (ah, SA[0.0, 0.0, 0.0]),
                (ts, r_z),
                (gh, r_z),
                (zp, r_z),
                (sp, r_z),
                (iw, r_z),
                (ps, SA[0.0, 0.0, 2.0e11]),
                (sph, SA[0.0, 0.0, 0.0]),
                (tok_c, r_z),
                (tok_p, SA[1.0, 0.0, 0.0]),
                (fad, r_mid),
            )
            JB = ForwardDiff.jacobian(field, pt)
            @test all(isfinite, JB)
            @test abs(tr(JB)) <= 1.0e-11 * max(1.0, opnorm(JB, Inf))
        end

        for (field, pt) in (
                (loop_ax, r_z),
                (loop_ax, r_mid),
                (sol, r_z),
                (sol_thick, r_z),
                (ts, r_z),
                (gh, r_z),
                (zp, r_z),
                (sp, r_z),
                (iw, r_z),
                (sph, SA[0.0, 0.0, 0.0]),
                (fad, r_mid),
            )
            JA = ForwardDiff.jacobian(r -> vector_potential(field, r), pt)
            @test all(isfinite, JA)
            @test isapprox(
                _ad_curl(r -> vector_potential(field, r), pt),
                field(pt);
                rtol = 1.0e-10,
                atol = 1.0e-12
            )
        end
    end
end

@testset "Composition, Grid, Geometry, and Reconstructions" begin
    let d1 = Dipole(SA[0.0, 0.0, 1.0], SA[0.5, 0.0, 0.0]),
            u1 = UniformField(SA[0.0, 0.0, 2.0e-6]),
            h1 = HarrisSheet(1.5e-5, 0.4),
            R = @SMatrix(
            [
                cos(0.7) -sin(0.7) 0.0;
                sin(0.7) cos(0.7) 0.0;
                0.0 0.0 1.0
            ]
        ),
            base_f = d1 + 2.0 * u1 - 0.5 * h1,
            comp = rotate_field(translate_field(base_f, SA[0.2, -0.4, 0.3]), R),
            pt = SVector(0.8, -0.3, 0.25)

        JB = ForwardDiff.jacobian(comp, pt)
        @test all(isfinite, JB)
        @test abs(tr(JB)) <= 1.0e-12 * max(1.0, opnorm(JB, Inf))
        @test isapprox(
            _ad_curl(r -> vector_potential(comp, r), pt), comp(pt); rtol = 1.0e-11
        )
        @test isapprox(
            _ad_curl(comp, pt),
            Magnetostatics.μ₀ * current_density(comp, pt);
            rtol = 1.0e-11
        )
    end

    let xs = range(-1.0, 1.0; length = 11),
            ys = range(-1.0, 1.0; length = 11),
            zs = range(-1.0, 1.0; length = 11),
            B = zeros(3, length(xs), length(ys), length(zs))
        for k in eachindex(zs), j in eachindex(ys), i in eachindex(xs)
            x, y, z = xs[i], ys[j], zs[k]
            B[1, i, j, k] = 2x - 3y + 0.5z
            B[2, i, j, k] = -x + 1.5y + 4z
            B[3, i, j, k] = 0.25x - y + 2z
        end
        J_exact = @SMatrix [2.0 -3.0 0.5; -1.0 1.5 4.0; 0.25 -1.0 2.0]
        pt = SVector(0.37, -0.42, 0.61)
        for interp in (:linear, :cubic)
            gf = GridMagneticField(B, xs, ys, zs; interpolation = interp)
            @test isapprox(ForwardDiff.jacobian(gf, pt), J_exact; rtol = 1.0e-12)
            @test all(isfinite, ForwardDiff.jacobian(r -> current_density(gf, r), pt))
        end
    end

    let zp = ZPinch(20.0, 0.4),
            pt = SVector(0.8, 0.2, 0.1)
        @test all(isfinite, ForwardDiff.jacobian(r -> b_unit(zp, r), pt))
        @test all(isfinite, ForwardDiff.gradient(r -> field_magnitude(zp, r), pt))
        @test all(isfinite, ForwardDiff.jacobian(r -> grad_b(zp, r), pt))
        @test all(isfinite, ForwardDiff.jacobian(r -> magnetic_curvature(zp, r), pt))
        @test all(isfinite, ForwardDiff.jacobian(r -> drift_velocities(zp, r).v_drift, pt))
    end

    let dom_ax = (0.8, 1.7, -0.8, 0.8),
            sol_eq = SolovevEquilibrium(2.0, 1.5, 1.5, 1.2),
            pts_rz = reduce(
            vcat, [
                reshape([R, Z], 1, 2)
                    for R in range(0.85, 1.65; length = 6)
                    for Z in range(-0.6, 0.6; length = 5)
            ]
        ),
            B_rz = reduce(
            vcat, [
                reshape(collect(sol_eq(SVector(pts_rz[i, 1], 0.0, pts_rz[i, 2]))), 1, 3)
                    for i in axes(pts_rz, 1)
            ]
        ),
            fit_ax = reconstruct_axisymmetric_field(
            pts_rz, B_rz; domain = dom_ax, length_scale = 0.6, regularization = 1.0e-12
        ),
            pt_ax = SVector(1.25, 0.15, 0.12)
        JB = ForwardDiff.jacobian(fit_ax, pt_ax)
        @test all(isfinite, JB)
        @test abs(tr(JB)) <= 1.0e-11 * max(1.0, opnorm(JB, Inf))
        @test isapprox(
            _ad_curl(fit_ax, pt_ax),
            Magnetostatics.μ₀ * current_density(fit_ax, pt_ax);
            rtol = 1.0e-10
        )
        @test all(isfinite, ForwardDiff.gradient(r -> poloidal_flux(fit_ax, r), pt_ax))
    end

    let box = (-0.5, 0.5, -0.5, 0.5, -0.5, 0.5),
            probes = reduce(
            vcat, [
                reshape([x, y, z], 1, 3)
                    for x in range(-0.4, 0.4; length = 4)
                    for y in range(-0.4, 0.4; length = 4)
                    for z in range(-0.4, 0.4; length = 4)
            ]
        ),
            truth = Dipole(SA[2.0e3, -1.0e3, 3.0e3], SA[0.1, -0.2, 1.8]) +
            0.5 * HarrisSheet(1.0, 0.4),
            B_obs = reduce(
            vcat,
            [reshape(collect(truth(SVector{3}(p))), 1, 3) for p in eachrow(probes)]
        ),
            pt = SVector(0.12, -0.09, 0.15)
        for (model, kernel) in (
                (:div_free, :wendland_c2),
                (:div_free, :gaussian),
                (:curl_free, :wendland_c2),
            )
            fit3d = reconstruct_field_3d(
                probes, B_obs;
                domain = box, model = model, kernel = kernel,
                length_scale = 0.55, regularization = 1.0e-8
            )
            JB = ForwardDiff.jacobian(fit3d, pt)
            @test all(isfinite, JB)
            @test abs(tr(JB)) <= 1.0e-11 * max(1.0, opnorm(JB, Inf))
            @test isapprox(
                _ad_curl(fit3d, pt),
                Magnetostatics.μ₀ * current_density(fit3d, pt);
                rtol = 1.0e-10,
                atol = 1.0e-12
            )
            if !isempty(fit3d.centers)
                c1 = SVector(fit3d.centers[1, 1], fit3d.centers[1, 2], fit3d.centers[1, 3])
                @test all(isfinite, ForwardDiff.jacobian(fit3d, c1))
            end
        end
    end

    let n_gsr = 25,
            times = collect(range(0.0, 2.0; length = n_gsr)),
            v = SVector(1.0, 0.0, 0.0),
            ez = SVector(0.0, 0.0, 1.0),
            B_obs = reduce(
            vcat, [
                reshape([0.2 * (t - 1.0), cos(0.4 * (t - 1.0)), 0.3], 1, 3)
                    for t in times
            ]
        ),
            p_obs = fill(0.1, n_gsr),
            ys = collect(range(-0.15, 0.15; length = 7)),
            gsr = reconstruct_grad_shafranov(
            times, Matrix(B_obs), p_obs;
            crossing_velocity = v, invariant_axis = ez, y = ys,
            spatial_degree = 6, pressure_degree = 2, axial_degree = 2
        ),
            pt_gsr = SVector(0.05, 0.03, 0.0)
        JB = ForwardDiff.jacobian(gsr, pt_gsr)
        @test all(isfinite, JB)
        @test abs(tr(JB)) <= 1.0e-11 * max(1.0, opnorm(JB, Inf))
        @test all(isfinite, ForwardDiff.gradient(r -> flux_function(gsr, r), pt_gsr))
        @test all(isfinite, ForwardDiff.gradient(r -> thermal_pressure(gsr, r), pt_gsr))
        @test isapprox(
            _ad_curl(gsr, pt_gsr),
            Magnetostatics.μ₀ * current_density(gsr, pt_gsr);
            rtol = 1.0e-10,
            atol = 1.0e-15
        )
    end
end

@testset "Solvers, Inputs, Float32, and Parameters" begin
    let wire = Wire([SA[0.0, 0.0, -2.0], SA[0.0, 0.0, 2.0]], 2.0),
            loop = CurrentLoop(0.8, 3.0, SA[0.0, 0.0, 0.2], SA[0.0, 0.0, 1.0]),
            ms_tup = MultiSource(wire, loop),
            ms_vec = MultiSource([loop]),
            bc_s = ConductingSphere(SA[0.0, 0.0, 0.0], 0.3, 6, 12),
            mesh = compute_surface_current(bc_s, wire),
            bs = BiotSavart(),
            bs_soft = BiotSavart(; core_radius = 0.05),
            bs_wall = BiotSavart(ConductingWall(3, -1.0)),
            bs_pbc = BiotSavart(PeriodicBoundary((2.0, 2.0, 0.0); n_images = (1, 1, 0))),
            bs_pre = BiotSavart(PrescribedBoundary(r -> SA[0.0, 0.0, 1.0e-3])),
            vp = VectorPotential(),
            vp_wall = VectorPotential(ConductingWall(3, -1.0)),
            pt = SVector(0.35, -0.25, 0.4)

        for (solver, src) in (
                (bs, wire),
                (bs_soft, wire),
                (bs, loop),
                (bs_soft, loop),
                (bs, mesh),
                (bs, ms_tup),
                (bs, ms_vec),
                (bs, [loop]),
                (bs_wall, wire),
                (bs_wall, loop),
                (bs_pbc, loop),
                (bs_pre, wire),
            )
            JB = ForwardDiff.jacobian(r -> solve(solver, src, r), pt)
            @test all(isfinite, JB)
            @test abs(tr(JB)) <= 1.0e-11 * max(1.0, opnorm(JB, Inf))
        end

        @test all(
            isfinite,
            ForwardDiff.jacobian(r -> solve(bs_soft, loop, r), SA[0.0, 0.0, 0.5])
        )

        for (vp_sol, bs_sol, src) in (
                (vp, bs, wire),
                (vp, bs, loop),
                (vp, bs, ms_tup),
                (vp, bs, ms_vec),
                (vp, bs, [loop]),
                (vp_wall, bs_wall, wire),
                (vp_wall, bs_wall, loop),
            )
            JA = ForwardDiff.jacobian(r -> solve(vp_sol, src, r), pt)
            @test all(isfinite, JA)
            @test isapprox(
                _ad_curl(r -> solve(vp_sol, src, r), pt),
                solve(bs_sol, src, pt);
                rtol = 1.0e-10
            )
        end

        dip = Dipole(SA[0.0, 0.0, 1.0])
        J_sv = ForwardDiff.jacobian(dip, pt)
        J_vec = ForwardDiff.jacobian(dip, collect(pt))
        J_tup = ForwardDiff.jacobian(r -> dip((r[1], r[2], r[3])), pt)
        J_xyz = ForwardDiff.jacobian(r -> dip(r[1], r[2], r[3]), pt)
        @test J_vec ≈ J_sv
        @test J_tup ≈ J_sv
        @test J_xyz ≈ J_sv

        B_ip = zeros(3)
        J_ip = ForwardDiff.jacobian((du, u) -> dip(du, u), B_ip, pt)
        @test J_ip ≈ J_sv

        dip32 = Dipole(SA[0.0f0, 0.0f0, 1.0f0])
        pt32 = SA[0.5f0, 0.2f0, -0.3f0]
        J32 = ForwardDiff.jacobian(dip32, pt32)
        @test eltype(J32) === Float32
        JA32 = ForwardDiff.jacobian(r -> vector_potential(dip32, r), pt32)
        @test eltype(JA32) === Float32

        r_check = SA[0.5, 0.2, 0.3]
        @test ForwardDiff.derivative(M -> Dipole(SA[0.0, 0.0, M])(r_check)[3], 8.0e22) ≈
            Dipole(SA[0.0, 0.0, 1.0])(r_check)[3]
        @test ForwardDiff.derivative(B0 -> HarrisSheet(B0, 1.0)(r_check)[1], 2.0) ≈
            tanh(r_check[3])
        @test ForwardDiff.derivative(I -> Solenoid(0.5, 1.0, I)(r_check)[3], 100.0) ≈
            Solenoid(0.5, 1.0, 1.0)(r_check)[3]
    end
end

end # module test_forwarddiff
