"""
    FFTSolver{BC <: AbstractBoundary, C} <: AbstractSolver

Solver that uses Fast Fourier Transforms to compute the magnetic field `B` from a current
distribution `J`.

The default constructor `FFTSolver()` uses `PeriodicBoundary`.  Use
`FFTSolver(OpenBoundary())` for free-space convolution without periodic wrap-around.
An optional grid-specific cache can be preallocated with
`FFTSolver(bc, (Nx, Ny, Nz), spacing; T = Float64)`.
"""
struct FFTSolver{BC <: AbstractBoundary, C} <: AbstractSolver
    bc::BC
    cache::C
end

function FFTSolver(bc::BC = PeriodicBoundary(SVector(0.0, 0.0, 0.0))) where {
        BC <: AbstractBoundary,
    }
    return FFTSolver{BC, Nothing}(bc, nothing)
end

struct _PeriodicFFTCache{T <: AbstractFloat, PF, PI}
    grid_size::NTuple{3, Int}
    spacing::NTuple{3, T}
    kx::Vector{T}
    ky::Vector{T}
    kz::Vector{T}
    plan_fwd::PF
    plan_inv::PI
    J_k::Array{Complex{T}, 4}
    B_k::Array{Complex{T}, 4}
    out::Array{T, 4}
end

struct _OpenFFTCache{T <: AbstractFloat, PF, PI}
    grid_size::NTuple{3, Int}
    spacing::NTuple{3, T}
    plan_fwd::PF
    plan_inv::PI
    J_pad::Array{T, 4}
    J_k::Array{Complex{T}, 4}
    B_k::Array{Complex{T}, 4}
    K_k::Array{Complex{T}, 4}
    G_k::Array{Complex{T}, 3}
    out::Array{T, 4}
end

@inline function _normalize_fft_spacing(::Type{T}, dx::Real) where {T <: AbstractFloat}
    v = T(dx)
    return (v, v, v)
end

@inline function _normalize_fft_spacing(
        ::Type{T}, spacing::NTuple{3, Real}
    ) where {T <: AbstractFloat}
    return (T(spacing[1]), T(spacing[2]), T(spacing[3]))
end

function _build_open_kernels(
        ::Type{T}, Nx::Int, Ny::Int, Nz::Int, dx::T, dy::T, dz::T
    ) where {T <: AbstractFloat}
    Px, Py, Pz = 2 * Nx, 2 * Ny, 2 * Nz
    K = zeros(T, 3, Px, Py, Pz)
    G = zeros(T, Px, Py, Pz)
    dV = dx * dy * dz
    coeff = T(μ₀) * dV / T(4π)
    self_G = T(μ₀) * T(0.1894) * cbrt(dV)^2

    @inbounds for k in 0:(Pz - 1)
        rz = (k < Nz ? k : k - Pz) * dz
        nyq_z = k == Nz
        rz2 = rz * rz
        for j in 0:(Py - 1)
            ry = (j < Ny ? j : j - Py) * dy
            nyq_yz = nyq_z || (j == Ny)
            ryz2 = ry * ry + rz2
            for i in 0:(Px - 1)
                if i == 0 && j == 0 && k == 0
                    G[1, 1, 1] = self_G
                    continue
                end
                rx = (i < Nx ? i : i - Px) * dx
                r2 = rx * rx + ryz2
                inv_r = inv(sqrt(r2))
                G[i + 1, j + 1, k + 1] = coeff * inv_r
                if !(nyq_yz || i == Nx)
                    fac = coeff * inv_r / r2
                    K[1, i + 1, j + 1, k + 1] = fac * rx
                    K[2, i + 1, j + 1, k + 1] = fac * ry
                    K[3, i + 1, j + 1, k + 1] = fac * rz
                end
            end
        end
    end

    return rfft(K, (2, 3, 4)), rfft(G)
end

function FFTSolver(
        bc::PeriodicBoundary, grid_size::NTuple{3, Int},
        spacing::Union{Real, NTuple{3, Real}}; T::Type{<:AbstractFloat} = Float64
    )
    Nx, Ny, Nz = grid_size
    dx, dy, dz = _normalize_fft_spacing(T, spacing)
    kx = T(2π) .* T.(rfftfreq(Nx, inv(dx)))
    ky = T(2π) .* T.(fftfreq(Ny, inv(dy)))
    kz = T(2π) .* T.(fftfreq(Nz, inv(dz)))
    sample_real = Array{T, 4}(undef, 3, Nx, Ny, Nz)
    plan_fwd = plan_rfft(sample_real, (2, 3, 4))
    nkx = div(Nx, 2) + 1
    sample_k = Array{Complex{T}, 4}(undef, 3, nkx, Ny, Nz)
    plan_inv = plan_irfft(sample_k, Nx, (2, 3, 4))
    J_k = similar(sample_k)
    B_k = similar(sample_k)
    out = similar(sample_real)
    cache = _PeriodicFFTCache(
        grid_size, (dx, dy, dz), kx, ky, kz, plan_fwd, plan_inv, J_k, B_k, out
    )
    return FFTSolver{PeriodicBoundary, typeof(cache)}(bc, cache)
end

function FFTSolver(
        bc::OpenBoundary, grid_size::NTuple{3, Int},
        spacing::Union{Real, NTuple{3, Real}}; T::Type{<:AbstractFloat} = Float64
    )
    Nx, Ny, Nz = grid_size
    dx, dy, dz = _normalize_fft_spacing(T, spacing)
    Px, Py, Pz = 2 * Nx, 2 * Ny, 2 * Nz
    K_k, G_k = _build_open_kernels(T, Nx, Ny, Nz, dx, dy, dz)
    J_pad = zeros(T, 3, Px, Py, Pz)
    plan_fwd = plan_rfft(J_pad, (2, 3, 4))
    J_k = similar(K_k)
    B_k = similar(K_k)
    plan_inv = plan_irfft(B_k, Px, (2, 3, 4))
    out = similar(J_pad)
    cache = _OpenFFTCache(
        grid_size, (dx, dy, dz), plan_fwd, plan_inv, J_pad, J_k, B_k, K_k, G_k, out
    )
    return FFTSolver{OpenBoundary, typeof(cache)}(bc, cache)
end

function _validate_fft_grid(J::AbstractArray{<:Real, 4}, xs, ys, zs)
    size(J, 1) == 3 || throw(ArgumentError("First dimension of J must be 3 (components)"))
    _, Nx, Ny, Nz = size(J)
    (length(xs) == Nx && length(ys) == Ny && length(zs) == Nz) ||
        throw(ArgumentError("Grid coordinate lengths must match J dimensions"))
    return (step(xs), step(ys), step(zs))
end

@inline function _fill_periodic_B_k!(B_k, J_k, kx, ky, kz, mu0::T) where {T}
    im_mu0 = Complex{T}(0, mu0)
    nkx, Ny, Nz = length(kx), length(ky), length(kz)
    @inbounds for k in 1:Nz
        kzv = kz[k]
        kz2 = kzv * kzv
        for j in 1:Ny
            kyv = ky[j]
            kyz2 = kyv * kyv + kz2
            for i in 1:nkx
                kxv = kx[i]
                k_sq = kxv * kxv + yz_sq(kyz2)
                if k_sq == zero(T)
                    B_k[1, i, j, k] = zero(Complex{T})
                    B_k[2, i, j, k] = zero(Complex{T})
                    B_k[3, i, j, k] = zero(Complex{T})
                else
                    fac = im_mu0 / k_sq
                    jx = J_k[1, i, j, k]
                    jy = J_k[2, i, j, k]
                    jz = J_k[3, i, j, k]
                    B_k[1, i, j, k] = fac * (kyv * jz - kzv * jy)
                    B_k[2, i, j, k] = fac * (kzv * jx - kxv * jz)
                    B_k[3, i, j, k] = fac * (kxv * jy - kyv * jx)
                end
            end
        end
    end
    return B_k
end

@inline yz_sq(x) = x

@inline function _fill_periodic_A_k!(A_k, J_k, kx, ky, kz, mu0::T) where {T}
    nkx, Ny, Nz = length(kx), length(ky), length(kz)
    @inbounds for k in 1:Nz
        kzv = kz[k]
        kz2 = kzv * kzv
        for j in 1:Ny
            kyv = ky[j]
            kyz2 = kyv * kyv + kz2
            for i in 1:nkx
                kxv = kx[i]
                k_sq = kxv * kxv + kyz2
                if k_sq == zero(T)
                    A_k[1, i, j, k] = zero(Complex{T})
                    A_k[2, i, j, k] = zero(Complex{T})
                    A_k[3, i, j, k] = zero(Complex{T})
                else
                    fac = mu0 / k_sq
                    jx = J_k[1, i, j, k]
                    jy = J_k[2, i, j, k]
                    jz = J_k[3, i, j, k]
                    k_dot_J_over_k2 = (kxv * jx + kyv * jy + kzv * jz) / k_sq
                    A_k[1, i, j, k] = fac * (jx - kxv * k_dot_J_over_k2)
                    A_k[2, i, j, k] = fac * (jy - kyv * k_dot_J_over_k2)
                    A_k[3, i, j, k] = fac * (jz - kzv * k_dot_J_over_k2)
                end
            end
        end
    end
    return A_k
end

@inline function _fill_open_B_k!(B_k, J_k, K_k)
    _, nkx, Py, Pz = size(J_k)
    @inbounds for k in 1:Pz, j in 1:Py, i in 1:nkx
        jx = J_k[1, i, j, k]
        jy = J_k[2, i, j, k]
        jz = J_k[3, i, j, k]
        kx = K_k[1, i, j, k]
        ky = K_k[2, i, j, k]
        kz = K_k[3, i, j, k]
        B_k[1, i, j, k] = jy * kz - jz * ky
        B_k[2, i, j, k] = jz * kx - jx * kz
        B_k[3, i, j, k] = jx * ky - jy * kx
    end
    return B_k
end

@inline function _fill_open_A_k!(A_k, J_k, G_k)
    _, nkx, Py, Pz = size(J_k)
    @inbounds for k in 1:Pz, j in 1:Py, i in 1:nkx
        gk = G_k[i, j, k]
        A_k[1, i, j, k] = J_k[1, i, j, k] * gk
        A_k[2, i, j, k] = J_k[2, i, j, k] * gk
        A_k[3, i, j, k] = J_k[3, i, j, k] * gk
    end
    return A_k
end

"""
    solve(solver::FFTSolver, J::AbstractArray{T, 4}, spacing) where T
    solve(solver::FFTSolver, J::AbstractArray{T, 4}, xs, ys, zs) where T

Compute the magnetic field `B` from a discrete current distribution `J` using FFT.
`spacing` can be an isotropic `dx::Real` or a 3-tuple `(dx, dy, dz)`.
"""
function solve(
        solver::FFTSolver, J::AbstractArray{T, 4}, dx::Real
    ) where {T <: Real}
    return solve(solver, J, (dx, dx, dx))
end

function solve(
        solver::FFTSolver, J::AbstractArray{<:Real, 4},
        xs::AbstractRange, ys::AbstractRange, zs::AbstractRange
    )
    spacing = _validate_fft_grid(J, xs, ys, zs)
    return solve(solver, J, spacing)
end

function _check_fft_cache(cache, grid_size, sp, ::Type{FT}) where {FT}
    (cache.grid_size == grid_size && isapprox.(cache.spacing, sp) |> all) ||
        throw(ArgumentError("FFTSolver cache does not match J size or grid spacing"))
    eltype(cache.out) === FT ||
        throw(ArgumentError("FFTSolver cache element type does not match float(eltype(J))"))
    return nothing
end

function solve(
        solver::FFTSolver{<:PeriodicBoundary, Nothing},
        J::AbstractArray{T, 4}, spacing::NTuple{3, Real}
    ) where {T <: Real}
    size(J, 1) == 3 || throw(ArgumentError("First dimension of J must be 3 (components)"))
    FT = float(T)
    _, Nx, Ny, Nz = size(J)
    dx, dy, dz = _normalize_fft_spacing(FT, spacing)
    kx = FT(2π) .* FT.(rfftfreq(Nx, inv(dx)))
    ky = FT(2π) .* FT.(fftfreq(Ny, inv(dy)))
    kz = FT(2π) .* FT.(fftfreq(Nz, inv(dz)))

    J_f = FT === T ? J : FT.(J)
    J_k = rfft(J_f, (2, 3, 4))
    B_k = similar(J_k)
    _fill_periodic_B_k!(B_k, J_k, kx, ky, kz, FT(μ₀))
    return irfft(B_k, Nx, (2, 3, 4))
end

function solve(
        solver::FFTSolver{<:PeriodicBoundary, C},
        J::AbstractArray{T, 4}, spacing::NTuple{3, Real}
    ) where {T <: Real, C <: _PeriodicFFTCache}
    size(J, 1) == 3 || throw(ArgumentError("First dimension of J must be 3 (components)"))
    FT = float(T)
    _, Nx, Ny, Nz = size(J)
    sp = _normalize_fft_spacing(FT, spacing)
    cache = solver.cache
    _check_fft_cache(cache, (Nx, Ny, Nz), sp, FT)

    J_f = FT === T ? J : FT.(J)
    mul!(cache.J_k, cache.plan_fwd, J_f)
    _fill_periodic_B_k!(cache.B_k, cache.J_k, cache.kx, cache.ky, cache.kz, FT(μ₀))
    mul!(cache.out, cache.plan_inv, cache.B_k)
    return copy(cache.out)
end

function solve(
        ::FFTSolver{OpenBoundary, Nothing},
        J::AbstractArray{T, 4}, spacing::NTuple{3, Real}
    ) where {T <: Real}
    size(J, 1) == 3 || throw(ArgumentError("First dimension of J must be 3 (components)"))
    FT = float(T)
    _, Nx, Ny, Nz = size(J)
    dx, dy, dz = _normalize_fft_spacing(FT, spacing)
    Px, Py, Pz = 2 * Nx, 2 * Ny, 2 * Nz

    K_k, _ = _build_open_kernels(FT, Nx, Ny, Nz, dx, dy, dz)
    J_pad = zeros(FT, 3, Px, Py, Pz)
    @inbounds J_pad[:, 1:Nx, 1:Ny, 1:Nz] .= J

    J_k = rfft(J_pad, (2, 3, 4))
    B_k = similar(J_k)
    _fill_open_B_k!(B_k, J_k, K_k)
    B_full = irfft(B_k, Px, (2, 3, 4))
    return B_full[:, 1:Nx, 1:Ny, 1:Nz]
end

function solve(
        solver::FFTSolver{OpenBoundary, C},
        J::AbstractArray{T, 4}, spacing::NTuple{3, Real}
    ) where {T <: Real, C <: _OpenFFTCache}
    size(J, 1) == 3 || throw(ArgumentError("First dimension of J must be 3 (components)"))
    FT = float(T)
    _, Nx, Ny, Nz = size(J)
    sp = _normalize_fft_spacing(FT, spacing)
    cache = solver.cache
    _check_fft_cache(cache, (Nx, Ny, Nz), sp, FT)

    fill!(cache.J_pad, zero(FT))
    @inbounds cache.J_pad[:, 1:Nx, 1:Ny, 1:Nz] .= J
    mul!(cache.J_k, cache.plan_fwd, cache.J_pad)
    _fill_open_B_k!(cache.B_k, cache.J_k, cache.K_k)
    mul!(cache.out, cache.plan_inv, cache.B_k)
    return cache.out[:, 1:Nx, 1:Ny, 1:Nz]
end

"""
    vector_potential(solver::FFTSolver, J::AbstractArray{T, 4}, spacing) where T
    vector_potential(solver::FFTSolver, J::AbstractArray{T, 4}, xs, ys, zs) where T

Compute the Coulomb-gauge magnetic vector potential `A` (`∇ ⋅ A = 0`) on a grid.
"""
function vector_potential(
        solver::FFTSolver, J::AbstractArray{T, 4}, dx::Real
    ) where {T <: Real}
    return vector_potential(solver, J, (dx, dx, dx))
end

function vector_potential(
        solver::FFTSolver, J::AbstractArray{<:Real, 4},
        xs::AbstractRange, ys::AbstractRange, zs::AbstractRange
    )
    spacing = _validate_fft_grid(J, xs, ys, zs)
    return vector_potential(solver, J, spacing)
end

function vector_potential(
        solver::FFTSolver{<:PeriodicBoundary, Nothing},
        J::AbstractArray{T, 4}, spacing::NTuple{3, Real}
    ) where {T <: Real}
    size(J, 1) == 3 || throw(ArgumentError("First dimension of J must be 3 (components)"))
    FT = float(T)
    _, Nx, Ny, Nz = size(J)
    dx, dy, dz = _normalize_fft_spacing(FT, spacing)
    kx = FT(2π) .* FT.(rfftfreq(Nx, inv(dx)))
    ky = FT(2π) .* FT.(fftfreq(Ny, inv(dy)))
    kz = FT(2π) .* FT.(fftfreq(Nz, inv(dz)))

    J_f = FT === T ? J : FT.(J)
    J_k = rfft(J_f, (2, 3, 4))
    A_k = similar(J_k)
    _fill_periodic_A_k!(A_k, J_k, kx, ky, kz, FT(μ₀))
    return irfft(A_k, Nx, (2, 3, 4))
end

function vector_potential(
        solver::FFTSolver{<:PeriodicBoundary, C},
        J::AbstractArray{T, 4}, spacing::NTuple{3, Real}
    ) where {T <: Real, C <: _PeriodicFFTCache}
    size(J, 1) == 3 || throw(ArgumentError("First dimension of J must be 3 (components)"))
    FT = float(T)
    _, Nx, Ny, Nz = size(J)
    sp = _normalize_fft_spacing(FT, spacing)
    cache = solver.cache
    _check_fft_cache(cache, (Nx, Ny, Nz), sp, FT)

    J_f = FT === T ? J : FT.(J)
    mul!(cache.J_k, cache.plan_fwd, J_f)
    _fill_periodic_A_k!(cache.B_k, cache.J_k, cache.kx, cache.ky, cache.kz, FT(μ₀))
    mul!(cache.out, cache.plan_inv, cache.B_k)
    return copy(cache.out)
end

function vector_potential(
        ::FFTSolver{OpenBoundary, Nothing},
        J::AbstractArray{T, 4}, spacing::NTuple{3, Real}
    ) where {T <: Real}
    size(J, 1) == 3 || throw(ArgumentError("First dimension of J must be 3 (components)"))
    FT = float(T)
    _, Nx, Ny, Nz = size(J)
    dx, dy, dz = _normalize_fft_spacing(FT, spacing)
    Px, Py, Pz = 2 * Nx, 2 * Ny, 2 * Nz

    _, G_k = _build_open_kernels(FT, Nx, Ny, Nz, dx, dy, dz)
    J_pad = zeros(FT, 3, Px, Py, Pz)
    @inbounds J_pad[:, 1:Nx, 1:Ny, 1:Nz] .= J

    J_k = rfft(J_pad, (2, 3, 4))
    A_k = similar(J_k)
    _fill_open_A_k!(A_k, J_k, G_k)
    A_full = irfft(A_k, Px, (2, 3, 4))
    return A_full[:, 1:Nx, 1:Ny, 1:Nz]
end

function vector_potential(
        solver::FFTSolver{OpenBoundary, C},
        J::AbstractArray{T, 4}, spacing::NTuple{3, Real}
    ) where {T <: Real, C <: _OpenFFTCache}
    size(J, 1) == 3 || throw(ArgumentError("First dimension of J must be 3 (components)"))
    FT = float(T)
    _, Nx, Ny, Nz = size(J)
    sp = _normalize_fft_spacing(FT, spacing)
    cache = solver.cache
    _check_fft_cache(cache, (Nx, Ny, Nz), sp, FT)

    fill!(cache.J_pad, zero(FT))
    @inbounds cache.J_pad[:, 1:Nx, 1:Ny, 1:Nz] .= J
    mul!(cache.J_k, cache.plan_fwd, cache.J_pad)
    _fill_open_A_k!(cache.B_k, cache.J_k, cache.G_k)
    mul!(cache.out, cache.plan_inv, cache.B_k)
    return cache.out[:, 1:Nx, 1:Ny, 1:Nz]
end
