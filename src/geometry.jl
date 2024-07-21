module Geometries

export Geometry

using ..Types

abstract type Geometry{T<:Real} end

# COV_EXCL_START

"Always true: this is an affine geometry."
function affine(g::Geometry{T})::Bool where {T<:Real} end

"Obtain the type of the reference element."
function type(g::Geometry{T})::GeometryType where {T<:Real} end

"Obtain number of corners of the corresponding reference element."
function corners(g::Geometry{T})::Integer where {T<:Real} end

"Obtain coordinates of the i-th corner."
function corner(g::Geometry{T}, i::Integer) where {T<:Real} end

"Obtain the centroid of the mapping's image."
function center(g::Geometry{T}) where {T<:Real} end

"""
  global(g,x)

Evaluate the local-to-global mapping.

# Arguments:
- `x::AbstractVector{S}`: local coordinate to map
"""
function toGlobal(g::Geometry{T}, x::AbstractVector{S}) where {T<:Real,S<:Real} end

"""
  local(g,X)

Evaluate the inverse mapping, i.e. the global-to-local mapping.

The returned local coordinate y minimizes
```
(global( y ) - x).two_norm()
```
on the entire affine hull of the reference element.  This degenerates
to the inverse map if the argument y is in the range of the map.

# Arguments:
- `X::AbstractVector{S}`: global coordinate to map
"""
function toLocal(g::Geometry{T}, X::AbstractVector{S}) where {T<:Real,S<:Real} end

raw"""
  integrationElement(g,x)

Obtain the integration element.

If the Jacobian of the mapping is denoted by ``J(x)``, the integration
integration element ``\mu(x)`` is given by
```math
\mu(x) = \sqrt{|\det (J^T(x) J(x))|}.
````

# Arguments:
- `x::AbstractVector{S}`: local coordinate to evaluate the integration element in
"""
function integrationElement(g::Geometry{T}, x::AbstractVector{S})::T where {T<:Real,S<:Real} end

"Obtain the volume of the element."
function volume(g::Geometry{T})::T where {T<:Real} end

"""
  jacobianTransposed(g,x)

Obtain the transposed of the Jacobian.

# Arguments:
- `x::AbstractVector{S}`: local coordinate to evaluate Jacobian in
"""
function jacobianTransposed(g::Geometry{T}, x::AbstractVector{S}) where {T<:Real,S<:Real} end

raw"""
  jacobianInverseTransposed(g,x)

Obtain the transposed of the Jacobian's inverse.

The Jacobian's inverse is defined as a pseudo-inverse. If we denote
the Jacobian by ``J(x)``, the following condition holds:
```math
J^{-1}(x) J(x) = I.
````

# Arguments:
- `x::AbstractVector{S}`: local coordinate to evaluate Jacobian in
"""
function jacobianInverseTransposed(g::Geometry{T}, x::AbstractVector{S}) where {T<:Real,S<:Real} end

"""
  jacobian(g,x)

Obtain the Jacobian.

# Arguments:
- `x::AbstractVector{S}`: local coordinate to evaluate Jacobian in
"""
function jacobian(g::Geometry{T}, x::AbstractVector{S}) where {T<:Real,S<:Real}
  transpose(jacobianTransposed(g,x))
end

raw"""
  jacobianInverse(g,x)

Obtain the Jacobian's inverse.

The Jacobian's inverse is defined as a pseudo-inverse. If we denote
the Jacobian by ``J(x)``, the following condition holds:
```math
J^{-1}(x) J(x) = I.
````

# Arguments:
- `x::AbstractVector{S}`: local coordinate to evaluate Jacobian in
"""
function jacobianInverse(g::Geometry{T}, x::AbstractVector{S}) where {T<:Real,S<:Real}
  transpose(jacobianInverseTransposed(g,x))
end

# COV_EXCL_STOP

end # module Geometries