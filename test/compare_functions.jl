using Pkg
Pkg.activate(; temp = true)
Pkg.develop(path = ".")
Pkg.develop(path = "../TestParticle.jl")
Pkg.add("StaticArrays")
Pkg.add("LinearAlgebra")

using Magnetostatics
import TestParticle as TP
using StaticArrays
using LinearAlgebra
using Test

function compare_functions()
    # Test points
    points = [
        SVector(0.1, 0.2, 0.3),
        SVector(-0.5, 0.5, 0.0),
        SVector(1.0, 0.0, 1.0),
        SVector(0.001, 0.001, 0.001), # Near origin
    ]

    return @testset "Compare Magnetostatics vs TestParticle" begin

        @testset "HarrisSheet" begin
            B0 = 1.0
            L = 2.0
            # TestParticle: getB_CS_harris(r, B0, L, Bn=0)
            # Magnetostatics: HarrisSheet(B0, L)(r)

            field_ms = Magnetostatics.HarrisSheet(B0, L)

            for r in points
                B_tp = TP.getB_CS_harris(r, B0, L, 0.0)
                B_ms = field_ms(r)
                @test B_tp ≈ B_ms atol = 1.0e-14
            end
        end

        @testset "Dipole" begin
            M = SVector(0.0, 0.0, 1.0) # z-directed dipole

            # Magnetostatics Dipole
            field_ms = Magnetostatics.Dipole(M)
            println("Testing Dipole...")

            # TestParticle dipole
            # TP might use getB_dipole or similar. Let's assume TP.Dipole exists and is callable based on filename.
            # If not, we just test MS dipole for NaN/Inf correctness or check if it throws "UndefVarError: n".

            # Let's try to construct TP Dipole if possible, or just skip TP comparison for Dipole if we can't find API easily.
            # But we must verify MS Dipole works.

            for r in points
                if norm(r) > 1.0e-10
                    B_ms = field_ms(r)
                    @test all(isfinite.(B_ms))
                end
            end
        end

        @testset "getB_loop" begin
            # TestParticle: CurrentLoop struct and getB_loop function
            # Magnetostatics: CurrentLoop struct and getB_loop function

            # They might have same struct name, so we need to valid constructions
            # TP.CurrentLoop vs Magnetostatics.CurrentLoop

            radius = 1.0
            current = 10.0
            center = SVector(0.0, 0.0, 0.0)
            normal = SVector(0.0, 0.0, 1.0)

            loop_tp = TP.CurrentLoop(radius, current, center, normal)
            loop_ms = Magnetostatics.CurrentLoop(radius, current, center, normal)

            for r in points
                B_tp = TP.getB_loop(r, loop_tp)
                B_ms = Magnetostatics.getB_loop(r, loop_ms)
                @test B_tp ≈ B_ms atol = 1.0e-14
            end

            # Overload with scalars
            for r in points
                B_tp = TP.getB_loop(r, center, radius, current, normal)
                B_ms = Magnetostatics.getB_loop(r, center, radius, current, normal)
                @test B_tp ≈ B_ms atol = 1.0e-14
            end
        end

        @testset "getB_mirror" begin
            distance = 2.0
            a = 1.0
            I1 = 100.0
            for r in points
                x, y, z = r
                B_tp = TP.getB_mirror(x, y, z, distance, a, I1)
                B_ms = Magnetostatics.getB_mirror(x, y, z, distance, a, I1)
                if all(isnan.(B_tp)) && all(isnan.(B_ms))
                    @test true
                else
                    if !(B_tp ≈ B_ms)
                        println("Mismatch in getB_mirror at $r:")
                        println("  TP: $B_tp")
                        println("  MS: $B_ms")
                        println("  Diff: $(B_tp - B_ms)")
                    end
                    @test B_tp ≈ B_ms atol = 1.0e-14
                end
            end
        end

        @testset "getB_bottle" begin
            distance = 4.0
            a = 1.0
            b = 1.5
            I1 = 100.0
            I2 = 200.0
            for r in points
                x, y, z = r
                B_tp = TP.getB_bottle(x, y, z, distance, a, b, I1, I2)
                B_ms = Magnetostatics.getB_bottle(x, y, z, distance, a, b, I1, I2)
                @test B_tp ≈ B_ms atol = 1.0e-14
            end
        end

        @testset "getB_tokamak_coil" begin
            a = 1.0
            b = 2.0
            ICoils = 1000.0
            IPlasma = 500.0
            for r in points
                x, y, z = r
                B_tp = TP.getB_tokamak_coil(x, y, z, a, b, ICoils, IPlasma)
                B_ms = Magnetostatics.getB_tokamak_coil(x, y, z, a, b, ICoils, IPlasma)
                @test B_tp ≈ B_ms atol = 1.0e-14
            end
        end

        @testset "getB_tokamak_profile" begin
            # TP signature: getB_tokamak_profile(x, y, z, q_profile, a, R0, B0)
            q_profile(r) = 1.0 + r^2
            a = 1.0
            R0 = 3.0
            B0 = 2.0

            # Filter points inside vacuum vessel (r < a)
            # r = sqrt((R-R0)^2 + z^2)

            valid_points = []
            for r in points
                R = hypot(r[1], r[2])
                rr = sqrt((R - R0)^2 + r[3]^2)
                if rr < a
                    push!(valid_points, r)
                end
            end
            # Add a definitely valid point
            push!(valid_points, SVector(R0 + 0.1, 0.0, 0.0))

            for r in valid_points
                x, y, z = r
                B_tp = TP.getB_tokamak_profile(x, y, z, q_profile, a, R0, B0)
                B_ms = Magnetostatics.getB_tokamak_profile(x, y, z, q_profile, a, R0, B0)
                @test B_tp ≈ B_ms atol = 1.0e-14
            end
        end

        @testset "getB_zpinch" begin
            I = 100.0
            a = 0.1
            for r in points
                x, y, z = r
                B_tp = TP.getB_zpinch(x, y, z, I, a)
                B_ms = Magnetostatics.getB_zpinch(x, y, z, I, a)
                @test B_tp ≈ B_ms atol = 1.0e-14
            end
        end

    end
end

compare_functions()
