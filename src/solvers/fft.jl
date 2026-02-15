using FFTW

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
    # fftfreq returns frequencies f = k / (2pi). We need k = 2pi * f.
    kx = 2π * fftfreq(Nx, 1 / dx)
    ky = 2π * fftfreq(Ny, 1 / dx)
    kz = 2π * fftfreq(Nz, 1 / dx)

    # Reshape for broadcasting
    kx_grid = reshape(kx, Nx, 1, 1)
    ky_grid = reshape(ky, 1, Ny, 1)
    kz_grid = reshape(kz, 1, 1, Nz)

    k_sq = zeros(eltype(kx), Nx, Ny, Nz)
    # Avoid loop for performance if possible, but simple loop is fine for setup
    @. k_sq = kx_grid^2 + ky_grid^2 + kz_grid^2

    # Handle singularity at k=0 (DC component)
    # If the total current is zero (sum(J)=0), then B(k=0) should be 0.
    # We set k^2 to non-zero (infinity) at index 1 to avoid division by zero,
    # effectively masking the DC component if we don't handle it explicitly.
    k_sq[1, 1, 1] = Inf

    # FFT of Current J
    # FFTW works on the last dimensions by default if we don't specify,
    # but here we have (3, Nx, Ny, Nz). We want to transform over dims 2,3,4.
    J_k = fft(J, (2, 3, 4))

    # Allocate B_k
    B_k = similar(J_k)

    # Calculation: B_k = i * mu0 * (k x J_k) / k^2
    im_mu0 = im * μ₀

    # It's more efficient to loop since we are dealing with vector components per grid point
    # Optimization: Use broadcasting or explicit loops

    # Let's perform the cross product and scaling
    # k x J = (ky Jz - kz Jy, kz Jx - kx Jz, kx Jy - ky Jx)

    for k in 1:Nz, j in 1:Ny, i in 1:Nx
        # Wave vectors at this grid point
        k_x = kx[i]
        k_y = ky[j]
        k_z = kz[k]

        kk = k_sq[i, j, k]

        # J_k at this grid point
        Jx = J_k[1, i, j, k]
        Jy = J_k[2, i, j, k]
        Jz = J_k[3, i, j, k]

        # Cross product k x J
        cross_x = k_y * Jz - k_z * Jy
        cross_y = k_z * Jx - k_x * Jz
        cross_z = k_x * Jy - k_y * Jx

        # Scale by i * mu0 / k^2
        factor = im_mu0 / kk

        B_k[1, i, j, k] = factor * cross_x
        B_k[2, i, j, k] = factor * cross_y
        B_k[3, i, j, k] = factor * cross_z
    end

    # Inverse FFT to get B in real space
    B = ifft(B_k, (2, 3, 4))

    # Return real part (imaginary part should be numerical noise for real inputs)
    return real.(B)
end
