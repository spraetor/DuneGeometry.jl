module AffineGeometries

export AffineGeometry

using LinearAlgebra: pinv, det
using ArgCheck

# Import other sub-modules
using ..Geometries
using ..ReferenceElements


struct AffineGeometry{T<:Real} <: Geometry{T}
  # public
  mydimension::Int
  coorddimension::Int

  # private
  refElement_::ReferenceElement{T}
  origin_::Vector{T}
  jacobianTransposed_::Array{T,2}
  jacobianInverseTransposed_::Array{T,2}
  integrationElement_::T
end

"Create affine geometry from reference element, one vertex, and the Jacobian matrix."
function AffineGeometry(refElement::ReferenceElement{T}, origin::AbstractVector{T},
                        jt::AbstractArray{T,2}) where {T<:Real}
  @argcheck size(jt,1) == length(origin)
  @argcheck size(jt,2) == refElement.dimension

  jit = LinearAlgebra.pinv(jt)
  integrationElement = abs(LinearAlgebra.det(jt))
  return AffineGeometry(refElement.dimension, length(origin), refElement, origin, jt, jit, integrationElement)
end

# "Create affine geometry from GeometryType, one vertex, and the Jacobian matrix."
# function AffineGeometry(gt::GeometryType, origin::AbstractVector{T},
#                         jt::AbstractArray{T,2}) where {T<:Real}
#   return AffineGeometry(ReferenceElement{T}(gt), origin, jt)
# end

"Create affine geometry from reference element and a vector of vertex coordinates."
function AffineGeometry(refElement::ReferenceElement{T}, coordVector::AbstractVector{C}) where {T<:Real, S<:Real, C<:AbstractVector{S}}
  origin = coordVector[1]
  jt = T[ coordVector[i+1] - origin for i = 1:refElement.dimension ]
  return AffineGeometry(refElement, origin, jt)
end

# "Create affine geometry from GeometryType and a vector of vertex coordinates."
# function AffineGeometry(gt::GeometryType, coordVector::AbstractVector{C}) where {T<:Real, C<:AbstractVector{T}}
#   return AffineGeometry(ReferenceElement{T}(gt), coordVector)
# end


"Always true: this is an affine geometry."
function affine(g::AffineGeometry{T})::Bool where {T<:Real}
  true
end

"Obtain the type of the reference element."
function type(g::AffineGeometry{T})::GeometryType where {T<:Real}
  type(g.refElement_)
end

"Obtain number of corners of the corresponding reference element."
function corners(g::AffineGeometry{T}) where {T<:Real}
  size(g.refElement_, g.mydimension)
end

"Obtain coordinates of the i-th corner."
function corner(g::AffineGeometry{T}, i::Integer) where {T<:Real}
  toGlobal(g, position(g.refElement_, i, mydimension))
end

"Obtain the centroid of the mapping's image."
function center(g::AffineGeometry{T}) where {T<:Real}
  toGlobal(position(g.refElement_, 0, 0))
end

"""
  toGlobal(g,x)

Evaluate the local-to-global mapping.

# Arguments:
- `x::AbstractVector{S}`: local coordinate to map
"""
function toGlobal(g::AffineGeometry{T}, x::AbstractVector{S}) where {T<:Real,S<:Real}
  g.origin_ + g.jacobianTransposed_' * x
end

"""
  toLocal(g,X)

Evaluate the inverse mapping, i.e. the global-to-local mapping.

The returned local coordinate y minimizes
```
(toGlobal( y ) - x).two_norm()
```
on the entire affine hull of the reference element.  This degenerates
to the inverse map if the argument y is in the range of the map.

# Arguments:
- `X::AbstractVector{S}`: global coordinate to map
"""
function toLocal(g::AffineGeometry{T}, X::AbstractVector{S}) where {T<:Real,S<:Real}
  g.jacobianInverseTransposed_' * (X - g.origin_)
end

raw"""
  integrationElement(g,x)

Obtain the integration element.

If the Jacobian of the mapping is denoted by ``J(x)``, the integration
integration element ``\mu(x)`` is given by

```math
\mu(x) = \sqrt{|\det (J^T(x) J(x))|}.
```

# Arguments:
- `x::AbstractVector{S}`: local coordinate to evaluate the integration element in
"""
function integrationElement(g::AffineGeometry{T}, x::AbstractVector{S}) where {T<:Real,S<:Real}
  g.integrationElement_;
end

"Obtain the volume of the element."
function volume(g::AffineGeometry{T}) where {T<:Real}
  g.integrationElement_ * volume(g.refElement_)
end

"""
  jacobianTransposed(g,x)

Obtain the transposed of the Jacobian.

# Arguments:
- `x::AbstractVector{S}`: local coordinate to evaluate Jacobian in
"""
function jacobianTransposed(g::AffineGeometry{T}, x::AbstractVector{S}) where {T<:Real,S<:Real}
  g.jacobianTransposed_
end

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
function jacobianInverseTransposed(g::AffineGeometry{T}, x::AbstractVector{S}) where {T<:Real,S<:Real}
  g.jacobianInverseTransposed_
end

"""
  jacobian(g,x)

Obtain the Jacobian.

# Arguments:
- `x::AbstractVector{S}`: local coordinate to evaluate Jacobian in
"""
function jacobian(g::AffineGeometry{T}, x::AbstractVector{S}) where {T<:Real,S<:Real}
  transpose(jacobianTransposed_)
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
function jacobianInverse(g::AffineGeometry{T}, x::AbstractVector{S}) where {T<:Real,S<:Real}
  transpose(jacobianInverseTransposed_)
end

"Obtain the reference element the geometry is defined on."
function referenceElement(g::AffineGeometry{T})::ReferenceElement{T} where {T<:Real}
  g.refElement_
end

end # module AffineGeometries
