@setup_workload begin
    let r = SVector(1.5, -0.4, 0.7),
            u = UniformField(SVector(0.0, 0.0, 1.0)),
            d = Dipole(SVector(0.0, 0.0, 1.0)),
            ld = LineDipole2D(SVector(0.0, -1.0); k = 1.0),
            lp = CurrentLoop(1.0, 1.0, SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 1.0)),
            cl = CurrentLoopAnalytic(lp; core_radius = 0.05),
            ws = WireSegment(SVector(-1.0, 0.0, 0.0), SVector(1.0, 0.0, 0.0), 1.0),
            tf = ToroidalField(1.0; core_radius = 0.05),
            g = [0.0 0.0; 1.0 0.0],
            h = zeros(2, 2),
            sh = SphericalHarmonicField(g, h, 1.0),
            sh_ws = SphericalHarmonicWorkspace(2)
        @compile_workload begin
            u(r); d(r); ld(r); cl(r); ws(r); tf(r)
            field_and_potential(d, r)
            field_and_potential(cl, r)
            field_and_potential(ws, r)
            field_and_potential(tf, r)
            evaluate!(sh_ws, sh, r)
            (u + d)(r)
        end
    end
end
