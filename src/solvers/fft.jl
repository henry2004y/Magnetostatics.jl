"""
    FFTSolver <: AbstractSolver

Solver that uses Fast Fourier Transforms to compute the magnetic field `B` from a current distribution `J`.
This method assumes periodic boundary conditions or sufficient padding.
"""
struct FFTSolver <: AbstractSolver end

"""
    solve(solver::FFTSolver, J::AbstractArray{T, 4}, dx::Real) where T

Compute the magnetic field `B` from a discrete current distribution `J` using FFT.

# Arguments
- `J`: 4D array of size (3, Nx, Ny, Nz) representing the current density components (Jx, Jy, Jz).
- `dx`: Grid spacing (assumed isotropic for now).

# Returns
- `B`: 4D array of size (3, Nx, Ny, Nz) representing the magnetic field components (Bx, By, Bz).
"""
function solve(::FFTSolver, J::AbstractArray{T, 4}, dx::Real) where {T}
    @assert size(J, 1) == 3 "First dimension of J must be 3 (components)"
    _, Nx, Ny, Nz = size(J)

    # k-vectors (angular wavenumbers)
    kx = 2π * fftfreq(Nx, 1 / dx)
    ky = 2π * fftfreq(Ny, 1 / dx)
    kz = 2π * fftfreq(Nz, 1 / dx)

    # FFT of Current J
    # FFTW works on the last dimensions by default, but here we have (3, Nx, Ny, Nz).
    J_k = fft(J, (2, 3, 4))

    # Allocate B_k
    B_k = similar(J_k)

    im_mu0 = im * μ₀

    @inbounds for k in 1:Nz, j in 1:Ny, i in 1:Nx
        # Wave vectors at this grid point
        k_vec = SVector(kx[i], ky[j], kz[k])
        k_sq = sum(abs2, k_vec)

        if k_sq == 0
            # Handle singularity at k=0 (DC component)
            B_k[:, i, j, k] .= 0
            continue
        end

        # J_k at this grid point
        J_vec = SVector(J_k[1, i, j, k], J_k[2, i, j, k], J_k[3, i, j, k])

        # Calculation: B_k = i * mu0 * (k x J_k) / k^2
        factor = im_mu0 / k_sq
        B_k[:, i, j, k] = factor * cross(k_vec, J_vec)
    end

    # Inverse FFT to get B in real space
    B = ifft(B_k, (2, 3, 4))

    # Return real part (imaginary part should be numerical noise for real inputs)
    return real.(B)
end
