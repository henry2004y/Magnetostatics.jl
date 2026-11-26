"""
    SphericalHarmonicWorkspace{T}

Reusable scratch buffers for zero-allocation evaluation of [`SphericalHarmonicField`](@ref).
Create one workspace per thread for concurrent evaluations.
"""
struct SphericalHarmonicWorkspace{T <: AbstractFloat}
    nmax::Int
    rnm::Matrix{T}
    aorn::Vector{T}
    xpcos::Vector{T}
    xpsin::Vector{T}
end

function SphericalHarmonicWorkspace{T}(nmax::Int) where {T <: AbstractFloat}
    nmax >= 1 || throw(ArgumentError("nmax must be at least 1"))
    return SphericalHarmonicWorkspace{T}(
        nmax,
        zeros(T, nmax + 1, nmax + 1),
        zeros(T, nmax + 2),
        zeros(T, nmax),
        zeros(T, nmax)
    )
end

function SphericalHarmonicWorkspace(nmax::Int; T::Type{<:AbstractFloat} = Float64)
    return SphericalHarmonicWorkspace{T}(nmax)
end

"""
    SphericalHarmonicField{T} <: AbstractMagneticField

Immutable, thread-safe Schmidt semi-normalized spherical-harmonic magnetic field
(BATSRUS / planetary crustal convention):
`V = r_ref ∑_{n=0}^{N} (r_ref / r_eff)^{n+1} ∑_{m=0}^n P_n^m(cosθ)(g_{nm} cos mφ + h_{nm} sin mφ)`,
`B = -∇V`, with `r_eff = max(r, r_ref)` clamped at the reference surface `r_ref`.

Use `evaluate!(ws, field, r)` with a preallocated [`SphericalHarmonicWorkspace`](@ref)
for allocation-free hot loops.
"""
struct SphericalHarmonicField{T <: AbstractFloat} <: AbstractMagneticField
    center::SVector{3, T}
    r_ref::T
    nmax::Int
    g::Matrix{T}
    h::Matrix{T}

    function SphericalHarmonicField(
            g::AbstractMatrix,
            h::AbstractMatrix,
            r_ref::Real;
            center::AbstractVector = SVector(0, 0, 0)
        )
        n = size(g, 1)
        n >= 1 || throw(ArgumentError("coefficient matrices must be non-empty"))
        size(g) == (n, n) == size(h) ||
            throw(ArgumentError("coefficient matrices must be square nmax × nmax"))
        r_ref > 0 || throw(ArgumentError("r_ref must be positive"))
        T = float(promote_type(eltype(g), eltype(h), typeof(r_ref), eltype(center)))
        return new{T}(
            SVector{3, T}(center[1], center[2], center[3]),
            T(r_ref),
            n,
            Matrix{T}(g),
            Matrix{T}(h)
        )
    end
end

"""
    evaluate!(ws::SphericalHarmonicWorkspace, field::SphericalHarmonicField, r) -> SVector{3}

Evaluate `field` at 3D point `r` using preallocated scratch buffers `ws`. Allocation-free
when called with arrays of matching size.
"""
function evaluate!(
        ws::SphericalHarmonicWorkspace{FT},
        field::SphericalHarmonicField{T},
        r
    ) where {FT <: AbstractFloat, T <: AbstractFloat}
    @boundscheck length(r) >= 3 ||
        throw(ArgumentError("Input must have at least 3 elements."))
    n = field.nmax
    ws.nmax >= n || throw(ArgumentError("workspace nmax $(ws.nmax) < field nmax $n"))
    NN = n - 1
    @inbounds begin
        xc = FT(r[1]) - FT(field.center[1])
        yc = FT(r[2]) - FT(field.center[2])
        zc = FT(r[3]) - FT(field.center[3])
    end
    r_mag = sqrt(xc * xc + yc * yc + zc * zc)
    r_mag <= zero(FT) && return zero(SVector{3, FT})

    reff = max(r_mag, FT(field.r_ref))
    ct = zc / r_mag
    st = sqrt(max(zero(FT), one(FT) - ct * ct))
    phi = atan(yc, xc)
    sinmin = FT(1.0e-6)

    xpcos = ws.xpcos
    xpsin = ws.xpsin
    aorn = ws.aorn
    R = ws.rnm
    g = field.g
    h = field.h
    @inbounds begin
        for im in 0:NN
            xpcos[im + 1] = cos(im * phi)
            xpsin[im + 1] = sin(im * phi)
        end
        ratio = FT(field.r_ref) / reff
        aorn[1] = one(FT)
        for k in 1:(NN + 2)
            aorn[k + 1] = ratio * aorn[k]
        end
        fill!(R, zero(FT))
        R[1, 1] = one(FT)
        if NN + 1 >= 1
            R[2, 1] = ct
        end
        for nn in 1:NN
            if nn == 1
                R[nn + 1, nn + 1] = st * R[nn, nn]
            else
                R[nn + 1, nn + 1] = sqrt(FT(nn - 0.5) / nn) * st * R[nn, nn]
            end
            if nn + 1 <= NN + 1
                R[nn + 2, nn + 1] = ct * sqrt(FT(2 * nn + 1)) * R[nn + 1, nn + 1]
            end
        end
        for m in 0:NN, l in (m + 2):NN
            R[l + 1, m + 1] = (
                ct * FT(2 * l - 1) * R[l, m + 1] -
                    R[l - 1, m + 1] * sqrt(FT((l + m - 1) * (l - m - 1)))
            ) / sqrt(FT(l * l - m * m))
        end

        br = zero(FT)
        bth = zero(FT)
        bph = zero(FT)
        for m in 0:NN, nn in m:NN
            rnm = R[nn + 1, m + 1]
            rnm1 = m + 1 <= NN + 1 ? R[nn + 1, m + 2] : zero(FT)
            dR = if m == 0
                -sqrt(FT((nn + 1) * nn) / 2) * rnm1
            elseif st <= sinmin
                -sqrt(FT((nn + m + 1) * (nn - m))) * rnm1
            else
                m * ct * rnm / st - sqrt(FT((nn + m + 1) * (nn - m))) * rnm1
            end
            cc = FT(g[nn + 1, m + 1])
            dd = FT(h[nn + 1, m + 1])
            cd = cc * xpcos[m + 1] + dd * xpsin[m + 1]
            a = aorn[nn + 3]
            br += (nn + 1) * a * rnm * cd
            bth -= a * dR * cd
            if st > sinmin
                bph -= a * rnm * m / st * (-cc * xpsin[m + 1] + dd * xpcos[m + 1])
            end
        end
        cp = n >= 2 ? xpcos[2] : cos(phi)
        sp = n >= 2 ? xpsin[2] : sin(phi)
        return SVector{3, FT}(
            br * (st * cp) + bth * (ct * cp) - bph * sp,
            br * (st * sp) + bth * (ct * sp) + bph * cp,
            br * ct - bth * st
        )
    end
end

function (field::SphericalHarmonicField{T})(r) where {T}
    FT = _promote_float(T, r)
    if FT <: AbstractFloat
        ws = SphericalHarmonicWorkspace{FT}(field.nmax)
        return evaluate!(ws, field, r)
    end
    # ForwardDiff Dual inputs fall back without allocating a float workspace.
    return _eval_sh_generic(field, r, FT)
end

function _eval_sh_generic(field::SphericalHarmonicField{T}, r, ::Type{FT}) where {T, FT}
    @boundscheck length(r) >= 3 ||
        throw(ArgumentError("Input must have at least 3 elements."))
    n = field.nmax
    NN = n - 1
    @inbounds begin
        xc = FT(r[1]) - FT(field.center[1])
        yc = FT(r[2]) - FT(field.center[2])
        zc = FT(r[3]) - FT(field.center[3])
    end
    r_mag = sqrt(xc * xc + yc * yc + zc * zc)
    r_mag <= zero(FT) && return zero(SVector{3, FT})
    reff = max(r_mag, FT(field.r_ref))
    ct = zc / r_mag
    st = sqrt(max(zero(FT), one(FT) - ct * ct))
    phi = atan(yc, xc)
    sinmin = FT(1.0e-6)

    R = zeros(FT, n + 1, n + 1)
    R[1, 1] = one(FT)
    if NN + 1 >= 1
        R[2, 1] = ct
    end
    for nn in 1:NN
        R[nn + 1, nn + 1] = nn == 1 ? st * R[nn, nn] :
            sqrt(FT(nn - 0.5) / nn) * st * R[nn, nn]
        if nn + 1 <= NN + 1
            R[nn + 2, nn + 1] = ct * sqrt(FT(2 * nn + 1)) * R[nn + 1, nn + 1]
        end
    end
    for m in 0:NN, l in (m + 2):NN
        R[l + 1, m + 1] = (
            ct * FT(2 * l - 1) * R[l, m + 1] -
                R[l - 1, m + 1] * sqrt(FT((l + m - 1) * (l - m - 1)))
        ) / sqrt(FT(l * l - m * m))
    end

    br = zero(FT)
    bth = zero(FT)
    bph = zero(FT)
    ratio = FT(field.r_ref) / reff
    @inbounds for m in 0:NN
        cm = cos(m * phi)
        sm = sin(m * phi)
        for nn in m:NN
            rnm = R[nn + 1, m + 1]
            rnm1 = m + 1 <= NN + 1 ? R[nn + 1, m + 2] : zero(FT)
            dR = if m == 0
                -sqrt(FT((nn + 1) * nn) / 2) * rnm1
            elseif st <= sinmin
                -sqrt(FT((nn + m + 1) * (nn - m))) * rnm1
            else
                m * ct * rnm / st - sqrt(FT((nn + m + 1) * (nn - m))) * rnm1
            end
            cc = FT(field.g[nn + 1, m + 1])
            dd = FT(field.h[nn + 1, m + 1])
            cd = cc * cm + dd * sm
            a = ratio^(nn + 2)
            br += (nn + 1) * a * rnm * cd
            bth -= a * dR * cd
            if st > sinmin
                bph -= a * rnm * m / st * (-cc * sm + dd * cm)
            end
        end
    end
    cp = cos(phi)
    sp = sin(phi)
    return SVector{3, FT}(
        br * (st * cp) + bth * (ct * cp) - bph * sp,
        br * (st * sp) + bth * (ct * sp) + bph * cp,
        br * ct - bth * st
    )
end

@inline function current_density(::SphericalHarmonicField{T}, r) where {T}
    return zero(SVector{3, _promote_float(T, r)})
end
