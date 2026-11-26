# Fused AGM evaluation of complete elliptic integrals K(m) and E(m) for 0 <= m < 1.
@inline _ellipke(m::Real) = _ellipke_impl(float(m))

@inline function _ellipke_impl(m::T) where {T <: Real}
    m < zero(T) && return (T(π) / 2, T(π) / 2)
    m >= one(T) && return (T(Inf), one(T))
    b = sqrt(one(T) - m)
    c = m / (2 * (one(T) + b))
    a = (one(T) + b) / 2
    b = sqrt(b)
    sum_val = muladd(c, c, m / 2)
    pow2 = T(2)
    for _ in 1:16
        abs(a - b) <= 4 * eps(T) * a && break
        c = (a - b) / 2
        sum_val += pow2 * c * c
        a_next = (a + b) / 2
        b = sqrt(a * b)
        a = a_next
        pow2 *= 2
    end
    K = T(π) / (a + b)
    E = K * (one(T) - sum_val)
    return K, E
end

# Bulirsch generalized complete elliptic integral cel(kc, p, c, s) without branch cuts.
@inline function _cel(kc::Real, p::Real, c::Real, s::Real)
    return _cel_impl(promote(float(kc), float(p), float(c), float(s))...)
end

@inline function _cel_impl(kc_in::T, p_in::T, c_in::T, s_in::T) where {T <: Real}
    kc = abs(kc_in)
    !(kc > zero(T)) && return T(NaN)
    e = kc
    m = one(T)
    p = p_in
    cc = c_in
    ss = s_in
    if p > zero(T)
        p = sqrt(p)
        ss = ss / p
    else
        f = kc * kc
        q = one(T) - f
        g = one(T) - p
        f = f - p
        q = (ss - cc * p) * q
        p = sqrt(f / g)
        cc = (cc - ss) / g
        ss = -q / (g * g * p) + cc * p
    end
    ca = sqrt(eps(T)) / 2
    for _ in 1:32
        f = cc
        cc = cc + ss / p
        g = e / p
        ss = 2 * (ss + f * g)
        p = g + p
        g = m
        m = kc + m
        abs(g - kc) <= g * ca && break
        kc = 2 * sqrt(e)
        e = kc * m
    end
    return T(π) / 2 * (ss + cc * m) / (m * (m + p))
end
