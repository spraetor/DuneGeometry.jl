export MultiLinearGeometry
export affine,type,corners,corner,center,toGlobal,toLocal,integrationElement,volume,jacobian,jacobianTransposed,jacobianInverse,jacobianInverseTransposed,referenceElement

using LinearAlgebra: pinv, det, SingularException
using ArgCheck

import ..MultiLinearGeometryImpl
import ..Utilities

struct MultiLinearGeometry{T<:Real} <: AbstractGeometry{T}
  # public
  mydimension::Int
  coorddimension::Int

  # private
  refElement_::ReferenceElement{T}
  corners_::Vector{Vector{T}}
end

"Create multilinear geometry from reference element and a vector of vertex coordinates."
function MultiLinearGeometry{T}(refElement::ReferenceElement{T}, corners::AbstractVector{C}) where {T<:Real, C<:AbstractVector{T}}
  @argcheck length(corners) > 0
  coorddim = length(corners[1])
  return MultiLinearGeometry{T}(refElement.dimension, coorddim, refElement, corners)
end

"Create multilinear geometry from GeometryType and a vector of vertex coordinates."
function MultiLinearGeometry{T}(gt::GeometryType, corners::AbstractVector{C}) where {T<:Real, C<:AbstractVector{T}}
  return MultiLinearGeometry{T}(ReferenceElement{T}(gt), corners)
end


"Is this mapping affine?"
function affine(g::MultiLinearGeometry{T})::Bool where {T<:Real}
  jt = zeros(T, g.mydimension, g.coorddimension)
  return MultiLinearGeometryImpl.affine(topologyId(g), g.mydimension, g.mydimension,
    MultiLinearGeometryImpl.CornerStorage(g.corners_), jt)
end

"Obtain the type of the reference element."
function type(g::MultiLinearGeometry{T})::GeometryType where {T<:Real}
  GeometryType(topologyId(g), g.mydimension)
end

"Obtain number of corners of the corresponding reference element."
function corners(g::MultiLinearGeometry{T}) where {T<:Real}
  length(g.corners_)
end

"Obtain coordinates of the i-th corner."
function corner(g::MultiLinearGeometry{T}, i::Integer) where {T<:Real}
  g.corners_[i]
end

"Obtain the centroid of the mapping's image."
function center(g::MultiLinearGeometry{T}) where {T<:Real}
  toGlobal(g, position(referenceElement(g), 1, 0))
end

"""
  toGlobal(g,x)

Evaluate the local-to-global mapping.

# Arguments:
- `x::AbstractVector{S}`: local coordinate to map
"""
function toGlobal(g::MultiLinearGeometry{T}, x::AbstractVector{S}) where {T<:Real,S<:Real}
  X = zeros(T, g.coorddimension)
  MultiLinearGeometryImpl.toGlobal(false, topologyId(g), g.mydimension, g.mydimension,
    MultiLinearGeometryImpl.CornerStorage(g.corners_), T(1), x, T(1), X)
  return X
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
function toLocal(g::MultiLinearGeometry{T}, X::AbstractVector{S}) where {T<:Real,S<:Real}
  tol = 16 * eps(T)
  x = position(referenceElement(g),1,0)
  dx = ones(eltype(x), length(x))
  affineMapping = affine(g)

  while sum(dx.^2) > tol
    # Newton's method: DF^n dx^n = F^n, x^{n+1} -= dx^n
    dglobal = toGlobal(g,x) - X
    try
      dx = jacobianTransposed(g,x) \ dglobal
    catch e
      if isa(e, SingularException)
        return ones(T, g.mydimension) .* typemax(T)
      else
        throw(e)
      end
    end

    # update x with correction
    x .-= dx
  end
  return x
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
function integrationElement(g::MultiLinearGeometry{T}, x::AbstractVector{S}) where {T<:Real,S<:Real}
  jt = jacobianTransposed(g,x)
  sqrt(abs(det(jt*jt')))
end

"Obtain the volume of the element."
function volume(g::MultiLinearGeometry{T}) where {T<:Real}
  # NOTE: This might be inexact for nonlinear maps, e.g. for surface elements
  integrationElement(g, position(referenceElement(g), 1, 0)) * volume(referenceElement(g));
end

"""
  jacobianTransposed(g,x)

Obtain the transposed of the Jacobian.

# Arguments:
- `x::AbstractVector{S}`: local coordinate to evaluate Jacobian in
"""
function jacobianTransposed(g::MultiLinearGeometry{T}, x::AbstractVector{S}) where {T<:Real,S<:Real}
  jt = zeros(T, g.mydimension, g.coorddimension)
  MultiLinearGeometryImpl.jacobianTransposed(false, topologyId(g), g.mydimension, g.mydimension,
    MultiLinearGeometryImpl.CornerStorage(g.corners_), T(1), x, T(1), jt)
  return jt
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
function jacobianInverseTransposed(g::MultiLinearGeometry{T}, x::AbstractVector{S}) where {T<:Real,S<:Real}
  Utilities.rightInvA(jacobianTransposed(g,x))
end

"""
  jacobian(g,x)

Obtain the Jacobian.

# Arguments:
- `x::AbstractVector{S}`: local coordinate to evaluate Jacobian in
"""
function jacobian(g::MultiLinearGeometry{T}, x::AbstractVector{S}) where {T<:Real,S<:Real}
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
function jacobianInverse(g::MultiLinearGeometry{T}, x::AbstractVector{S}) where {T<:Real,S<:Real}
  transpose(jacobianInverseTransposed(g,x))
end

"Obtain the reference element the geometry is defined on."
function referenceElement(g::MultiLinearGeometry{T})::ReferenceElement{T} where {T<:Real}
  g.refElement_
end


function topologyId(g::MultiLinearGeometry{T})::UInt32 where {T<:Real}
    type(g.refElement_).topologyId
end


