"""
    FFTSolver{BC <: AbstractBoundary} <: AbstractSolver

Solver that uses Fast Fourier Transforms to compute the magnetic field `B` from a current
distribution `J`.

The default constructor `FFTSolver()` uses `PeriodicBoundary`, which is inherent to the
FFT method.  Use `FFTSolver(OpenBoundary())` to enable zero-padding for non-periodic
(free-space) fields.
"""
struct FFTSolver{BC <: AbstractBoundary} <: AbstractSolver
    bc::BC
end

FFTSolver() = FFTSolver(PeriodicBoundary(SVector(0.0, 0.0, 0.0)))

"""
    solve(solver::FFTSolver, J::AbstractArray{T, 4}, dx::Real) where T

Compute the magnetic field `B` from a discrete current distribution `J` using FFT.

# Arguments
- `J`: 4D array of size (3, Nx, Ny, Nz) representing the current density components (Jx, Jy, Jz).
- `dx`: Grid spacing (assumed isotropic for now).

# Returns
- `B`: 4D array of size (3, Nx, Ny, Nz) representing the magnetic field components (Bx, By, Bz).
"""
function solve(
        ::FFTSolver{<:PeriodicBoundary}, J::AbstractArray{T, 4},
        dx::Real
    ) where {T}
    @assert size(J, 1) == 3 "First dimension of J must be 3 (components)"
    _, Nx, Ny, Nz = size(J)
    invdx = inv(dx)
    kx = 2π * fftfreq(Nx, invdx)
    ky = 2π * fftfreq(Ny, invdx)
    kz = 2π * fftfreq(Nz, invdx)

    J_k = fft(J, (2, 3, 4))
    B_k = similar(J_k)
    im_mu0 = im * μ₀

    @inbounds for k in 1:Nz, j in 1:Ny, i in 1:Nx
        k_vec = SVector(kx[i], ky[j], kz[k])
        k_sq = sum(abs2, k_vec)
        if k_sq == 0
            B_k[:, i, j, k] .= 0
            continue
        end
        J_vec = SVector(J_k[1, i, j, k], J_k[2, i, j, k], J_k[3, i, j, k])
        factor = im_mu0 / k_sq
        B_k[:, i, j, k] = factor * cross(k_vec, J_vec)
    end

    return real.(ifft(B_k, (2, 3, 4)))
end

"""
    solve(solver::FFTSolver{OpenBoundary}, J, dx)

Open-boundary (free-space) solver using Hockney's discrete Green's function
convolution method. The current density `J` is zero-padded to twice the size
in each dimension and convolved with the isolated 3D Green's function
`G(r) = 1/(4π r)`, eliminating periodic wrap-around artifacts.
"""
function solve(::FFTSolver{OpenBoundary}, J::AbstractArray{T, 4}, dx::Real) where {T}
    @assert size(J, 1) == 3 "First dimension of J must be 3 (components)"
    _, Nx, Ny, Nz = size(J)
    Px, Py, Pz = 2Nx, 2Ny, 2Nz
    invdx = inv(dx)

    G = zeros(T, Px, Py, Pz)
    inv_4pi = inv(T(4π))
    self_val = T(0.1894) / dx

    @inbounds for k in 0:(Pz - 1), j in 0:(Py - 1), i in 0:(Px - 1)
        rx = (i <= Nx ? i : i - Px) * dx
        ry = (j <= Ny ? j : j - Py) * dx
        rz = (k <= Nz ? k : k - Pz) * dx
        r_sq = rx^2 + ry^2 + rz^2
        if r_sq == 0
            G[i + 1, j + 1, k + 1] = self_val
        else
            G[i + 1, j + 1, k + 1] = inv_4pi / sqrt(r_sq)
        end
    end

    G_k = fft(G) .* (dx^3)

    J_pad = zeros(complex(T), 3, Px, Py, Pz)
    J_pad[:, 1:Nx, 1:Ny, 1:Nz] .= J

    kx = 2π * fftfreq(Px, invdx)
    ky = 2π * fftfreq(Py, invdx)
    kz = 2π * fftfreq(Pz, invdx)

    J_k = fft(J_pad, (2, 3, 4))
    B_k = similar(J_k)
    im_mu0 = im * T(μ₀)

    @inbounds for k in 1:Pz, j in 1:Py, i in 1:Px
        k_vec = SVector(kx[i], ky[j], kz[k])
        factor = im_mu0 * G_k[i, j, k]
        J_vec = SVector(J_k[1, i, j, k], J_k[2, i, j, k], J_k[3, i, j, k])
        B_k[:, i, j, k] = factor * cross(k_vec, J_vec)
    end

    B_full = real.(ifft(B_k, (2, 3, 4)))
    return B_full[:, 1:Nx, 1:Ny, 1:Nz]
end

