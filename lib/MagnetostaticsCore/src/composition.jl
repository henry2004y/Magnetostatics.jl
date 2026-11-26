"""
    CompositeField{Fs <: Tuple} <: AbstractMagneticField

Superposition of multiple `AbstractMagneticField` models evaluated by summing their
contributions at each query position. Prefer constructing via `+` and `-`.
"""
struct CompositeField{Fs <: Tuple} <: AbstractMagneticField
    fields::Fs

    function CompositeField{Fs}(fields::Fs) where {Fs <: Tuple}
        !isempty(fields) ||
            throw(ArgumentError("CompositeField requires at least one field"))
        return new{Fs}(fields)
    end
end

CompositeField(fields::Fs) where {Fs <: Tuple} = CompositeField{Fs}(fields)
CompositeField(fields::AbstractMagneticField...) = CompositeField(fields)

@inline _eval_composite_field(fields::Tuple{F}, r) where {F} = fields[1](r)
@inline function _eval_composite_field(fields::Tuple, r)
    return fields[1](r) + _eval_composite_field(Base.tail(fields), r)
end

@inline (field::CompositeField)(r) = _eval_composite_field(field.fields, r)

@inline _eval_composite_j(fields::Tuple{F}, r) where {F} = current_density(fields[1], r)
@inline function _eval_composite_j(fields::Tuple, r)
    return current_density(fields[1], r) + _eval_composite_j(Base.tail(fields), r)
end

@inline function current_density(field::CompositeField, r)
    return _eval_composite_j(field.fields, r)
end

@inline _eval_composite_a(fields::Tuple{F}, r) where {F} = vector_potential(fields[1], r)
@inline function _eval_composite_a(fields::Tuple, r)
    return vector_potential(fields[1], r) + _eval_composite_a(Base.tail(fields), r)
end

@inline function vector_potential(field::CompositeField, r)
    return _eval_composite_a(field.fields, r)
end

@inline function _eval_composite_ba(fields::Tuple{F}, r) where {F}
    return field_and_potential(fields[1], r)
end
@inline function _eval_composite_ba(fields::Tuple, r)
    B1, A1 = field_and_potential(fields[1], r)
    B2, A2 = _eval_composite_ba(Base.tail(fields), r)
    return (B1 + B2, A1 + A2)
end

@inline function field_and_potential(field::CompositeField, r)
    return _eval_composite_ba(field.fields, r)
end

@inline _eval_composite_psi(fields::Tuple{F}, r) where {F} = poloidal_flux(fields[1], r)
@inline function _eval_composite_psi(fields::Tuple, r)
    return poloidal_flux(fields[1], r) + _eval_composite_psi(Base.tail(fields), r)
end

@inline function poloidal_flux(field::CompositeField, r)
    return _eval_composite_psi(field.fields, r)
end

Base.:+(::NullField, ::NullField) = NullField()
Base.:+(::NullField, f::AbstractMagneticField) = f
Base.:+(f::AbstractMagneticField, ::NullField) = f
Base.:+(::NullField, f::CompositeField) = f
Base.:+(f::CompositeField, ::NullField) = f
Base.:+(f1::AbstractMagneticField, f2::AbstractMagneticField) = CompositeField((f1, f2))
function Base.:+(f1::CompositeField, f2::AbstractMagneticField)
    return CompositeField((f1.fields..., f2))
end
function Base.:+(f1::AbstractMagneticField, f2::CompositeField)
    return CompositeField((f1, f2.fields...))
end
function Base.:+(f1::CompositeField, f2::CompositeField)
    return CompositeField((f1.fields..., f2.fields...))
end

Base.:-(f1::AbstractMagneticField, f2::AbstractMagneticField) = f1 + (-f2)

"""
    ScaledField{T <: Real, F <: AbstractMagneticField} <: AbstractMagneticField

Scalar multiple `scale * field`. Prefer constructing via `*`, `/`, and unary `-`.
"""
struct ScaledField{T <: Real, F <: AbstractMagneticField} <: AbstractMagneticField
    scale::T
    field::F
end

@inline (f::ScaledField)(r) = f.scale * f.field(r)
@inline current_density(f::ScaledField, r) = f.scale * current_density(f.field, r)
@inline vector_potential(f::ScaledField, r) = f.scale * vector_potential(f.field, r)
@inline function field_and_potential(f::ScaledField, r)
    B, A = field_and_potential(f.field, r)
    return (f.scale * B, f.scale * A)
end
@inline poloidal_flux(f::ScaledField, r) = f.scale * poloidal_flux(f.field, r)

Base.:*(::Real, ::NullField) = NullField()
Base.:*(::NullField, ::Real) = NullField()
Base.:*(c::Real, f::ScaledField) = ScaledField(c * f.scale, f.field)
Base.:*(f::ScaledField, c::Real) = ScaledField(f.scale * c, f.field)
Base.:*(c::Real, f::AbstractMagneticField) = ScaledField(c, f)
Base.:*(f::AbstractMagneticField, c::Real) = ScaledField(c, f)
Base.:-(f::NullField) = f
Base.:-(f::ScaledField) = ScaledField(-f.scale, f.field)
Base.:-(f::AbstractMagneticField) = ScaledField(-1, f)
Base.:/(f::AbstractMagneticField, c::Real) = inv(c) * f

"""
    TransformedField{T <: Real, F <: AbstractMagneticField}

Rigidly translated and rotated magnetic field: evaluates `field` in local coordinates
`r_loc = rotation' * (r - origin)` and rotates output vectors back via `rotation`.
"""
struct TransformedField{T <: Real, F <: AbstractMagneticField} <: AbstractMagneticField
    field::F
    origin::SVector{3, T}
    rotation::SMatrix{3, 3, T, 9}

    function TransformedField(
            field::F, origin::AbstractVector,
            rotation::AbstractMatrix = SMatrix{3, 3, Bool, 9}(I)
        ) where {F <: AbstractMagneticField}
        T = float(promote_type(eltype(origin), eltype(rotation)))
        R = SMatrix{3, 3, T, 9}(rotation)
        orig = SVector{3, T}(origin[1], origin[2], origin[3])
        isapprox(R' * R, I; atol = 1.0e-6) ||
            throw(ArgumentError("rotation matrix must be orthogonal"))
        det(R) > zero(T) ||
            throw(ArgumentError("rotation matrix must be proper (det ≈ +1)"))
        return new{T, F}(field, orig, R)
    end
end

@inline function _local_point(f::TransformedField{T}, r) where {T}
    @boundscheck length(r) >= 3 ||
        throw(ArgumentError("Input must have at least 3 elements."))
    FT = _promote_float(T, r)
    @inbounds dr = SVector{3, FT}(
        r[1] - f.origin[1], r[2] - f.origin[2], r[3] - f.origin[3]
    )
    return f.rotation' * dr
end

@inline (f::TransformedField)(r) = f.rotation * f.field(_local_point(f, r))
@inline function current_density(f::TransformedField, r)
    return f.rotation * current_density(f.field, _local_point(f, r))
end
@inline function vector_potential(f::TransformedField, r)
    return f.rotation * vector_potential(f.field, _local_point(f, r))
end
@inline function field_and_potential(f::TransformedField, r)
    B, A = field_and_potential(f.field, _local_point(f, r))
    return (f.rotation * B, f.rotation * A)
end

"""
    translate_field(field, offset) -> AbstractMagneticField

Return a field displaced by `offset` so that evaluating the result at `r` evaluates `field`
at `r - offset`.
"""
function translate_field(field::AbstractMagneticField, offset::AbstractVector)
    T = float(eltype(offset))
    return TransformedField(field, SVector{3, T}(offset), SMatrix{3, 3, T, 9}(I))
end

translate_field(f::NullField, ::AbstractVector) = f

function translate_field(f::TransformedField{T}, offset::AbstractVector) where {T}
    FT = float(promote_type(T, eltype(offset)))
    new_orig = SVector{3, FT}(f.origin) + SVector{3, FT}(offset)
    return TransformedField(f.field, new_orig, SMatrix{3, 3, FT, 9}(f.rotation))
end

"""
    rotate_field(field, rotation; center = SA[0, 0, 0]) -> AbstractMagneticField

Rotate `field` by proper orthogonal matrix `rotation` about `center`.
"""
function rotate_field(
        field::AbstractMagneticField, rotation::AbstractMatrix;
        center::AbstractVector = SVector(0, 0, 0)
    )
    T = float(promote_type(eltype(rotation), eltype(center)))
    Q = SMatrix{3, 3, T, 9}(rotation)
    c = SVector{3, T}(center[1], center[2], center[3])
    origin = c - Q * c
    return TransformedField(field, origin, Q)
end

rotate_field(f::NullField, ::AbstractMatrix; center = SVector(0, 0, 0)) = f

function rotate_field(
        f::TransformedField{T}, rotation::AbstractMatrix;
        center::AbstractVector = SVector(0, 0, 0)
    ) where {T}
    FT = float(promote_type(T, eltype(rotation), eltype(center)))
    Q = SMatrix{3, 3, FT, 9}(rotation)
    c = SVector{3, FT}(center[1], center[2], center[3])
    new_orig = c + Q * (SVector{3, FT}(f.origin) - c)
    return TransformedField(f.field, new_orig, Q * SMatrix{3, 3, FT, 9}(f.rotation))
end
