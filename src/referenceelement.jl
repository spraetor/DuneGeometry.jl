# Some type aliases
Coordinate{T<:Real} = Vector{T}

struct ReferenceElementGeometry{T<:Real}

end


"""
    ReferenceElement{T}

This class provides access to geometric and topological properties of a
reference element.

This includes the number of subentities, the volume, and a method for checking
if a point is contained in the reference element.
The embedding of each subentity into the reference element is also
provided.
"""
struct ReferenceElement{T<:Real}
    dimension::Int

    "The reference element volume."
    volume_::T

    baryCenters_::Vector{Vector{Coordinate{T}}}
    integrationNormals_::Vector{Coordinate{T}}

    geometries_::Vector{Vector{ReferenceElementGeometry{T}}}

    info_::Vector{Vector{SubEntityInfo}}

    "Constructor"
    ReferenceElement(dim::Int) = new(dim, 0::T,                     # dimension, volume_
        Vector{Vector{Coordinate{T}}}(undef, dim+1),                # baryCenters_
        undef,                                                      # integrationNormals_
        Vector{Vector{ReferenceElementGeometry{T}}}(undef, dim+1),  # geometries_
        Vector{Vector{SubEntityInfo}}(undef, dim+1))                # info_
end


"""
    size(r, c)

Number of subentities of codimension `c` in the ReferenceElement `r`.

# Arguments
- `r::ReferenceElement`: The reference element.
- `c::Int`: Codimension whose size is requested.
"""
function size{T<:Real}(self::ReferenceElement{T}, c::Int)::Int
    @argcheck 0 <= c <= self.dimension
    return _size(self.info_[c+1])
end


"""
    size(r, i, c, cc)

Number of subentities of codimension `cc` of subentity `(i,c)`.

Denote by `E` the `i`-th subentity of codimension `c` of the given
reference element `r`. This method returns the number of subentities
of codimension `cc` of the reference element, that are also
a subentity of `E`. If `cc<c` this number is zero.

# Arguments
- `r::ReferenceElement`: The reference element.
- `i:Int`: The number of subentity `E` (`0 < i <= size(r,c)`)
- `c::Int`: Codimension of subentity `E` (`0 <= c <= dim(r)`)
- `cc::Int`: Codimension whose size is desired (`0 <= cc <= dim(r)`)
"""
function size{T<:Real}(self::ReferenceElement{T}, i::Int, c::Int, cc::Int)::Int
    @argcheck 0 < i <= size(self, c)
    return _size(self.info_[c+1][i], cc)
end


"""
    subEntity(r, i, c, ii, cc)

Obtain number of `ii`-th subentity with codim `cc` of `(i,c)`.

Denote by E the i-th subentity of codimension c of the current
reference element. And denote by S the ii-th subentity of codimension
(cc-c) of E. Then, S is a also a subentity of codimension cc of the current
reference element. This method returns the number of S with respect
to the current reference element.

# Arguments
- `r::ReferenceElement`: The reference element.
- `i::Int`: Number of subentity E (0 < i <= size( c ))
- `c::Int`: Codimension of subentity E
- `ii:Int`: Number of subentity S (with respect to E)
- `cc:Int`: Codimension of subentity S (c <= cc <= dim)
"""
function subEntity{T<:Real}(self::ReferenceElement{T}, i::Int, c::Int, ii::Int, cc::Int)::Int
    @argcheck 0 < i <= size(self, c)
    return _number(self.info_[c+1][i], ii, cc )
end


"""
    subEntities(r, i, c, cc)

Obtain the range of numbers of subentities with codim `cc` of `(i,c)`.

Denote by E the `i`-th subentity of codimension `c` of the current
reference element. This method returns a range of numbers of
all subentities of E with codimension `cc`. Notice that the sub-subentity
codimension as well as the numbers in the returned range are
given with respect to the reference element itself and not with
respect to E. For `0<=cc<c` this will return an empty range.

# Arguments
- `r::ReferenceElement`: The reference element.
- `i::Int`: Number of subentity E (0 < i <= size( c ))
- `c::Int`: Codimension of subentity E
- `cc:Int`: Codimension of subentity S (c <= cc <= dim)
"""
function subEntities{T<:Real}(self::ReferenceElement{T}, i::Int, c::Int, cc::Int)
    @argcheck 0 < i <= size(self, c)
    return _numbers(self.info_[c+1][i], cc)
end


"""
    type(r, i, c)

Obtain the type of subentity `(i,c)`.

Denote by E the `i`-th subentity of codimension `c` of the current
reference element. This method returns the GeometryType of E.

# Arguments
- `r::ReferenceElement`: The reference element.
- `i::Int`: Number of subentity E (0 < i <= size( c ))
- `c::Int`: Codimension of subentity E
"""
function type{T<:Real}(self::ReferenceElement{T}, i::Int, c::Int)::GeometryType
    @argcheck 0 < i <= size(self, c)
    return _type(self.info_[c+1][i])
end


"Obtain the type of this reference element"
type{T<:Real}(self::ReferenceElement{T})::GeometryType = type(self, 0, 0)


"""
    position(r, i, c)

Position of the barycenter of entity `(i,c)`.

Denote by E the `i`-th subentity of codimension `c` of the current
reference element. This method returns the coordinates of
the center of gravity of E within the current reference element.

# Arguments
- `r::ReferenceElement`: The reference element.
- `i::Int`: Number of subentity E (0 < i <= size( c ))
- `c::Int`: Codimension of subentity E
"""
function position{T<:Real}(self::ReferenceElement{T}, i::Int, c::Int)::Vector{T}
    @argcheck 0 <= c <= dim
    return self.baryCenters_[c+1][i]
end


"""
    checkInside(r, local)

Check if a coordinate is in the reference element.

This method returns true if the given local coordinate is within this
reference element.

# Arguments
- `x::DenseVector`: Coordinates of the point.
"""
function checkInside{T<:Real}(self::ReferenceElement{T}, x::DenseVector{T})::Bool
    return _checkInside(type(self).topologyId, dimension, x, 64 * eps(T))
end


"""
    geometry(r, i)

Obtain the embedding of subentity `(i,codim)` into the reference element

Denote by E the i-th subentity of codimension codim of the current
reference element. This method returns a \ref Dune::AffineGeometry
that maps the reference element of E into the current reference element.

# Arguments
- `i::Int`: number of subentity E (0 < i <= size( c ))
- `c::Int`: Codimension of subentity E
"""
function geometry{T<:Real}(self::ReferenceElement{T}, i::Int, c::Int)::ReferenceElementGeometry{T}
    return self.geometries_[c+1][i]
end


"Obtain the volume of the reference element."
volume{T<:Real}(self::ReferenceElement{T})::T = self.volume_


"""
    integrationOuterNormal(r, face)

Obtain the integration outer normal of the reference element.

The integration outer normal is the outer normal whose length coincides
with the face's integration element.

# Arguments
- `face::Int`: Index of the face, whose normal is desired.
"""
function integrationOuterNormal{T<:Real}(self::ReferenceElement{T}, face::Int)::Vector{T}
    @argcheck 0 < face <= Base.length(self.integrationNormals_)
    return self.integrationNormals_[face]
end


function initialize{T<:Real}(self::ReferenceElement{T}, topologyId::UInt)
  @argcheck topologyId < numTopologies(self.dim)

  # set up subentities
  dim = self.dimension
  for codim = 0:dim
    s = _size(topologyId, dim, codim)
    resize!(self.info_[codim+1], s)

    for i = 1:s
      initialize(self.info_[codim+1][i], topologyId, codim, i)
    end
  end

  # compute corners
  numVertices = size(self,dim)
  resize!(self.baryCenters_[dim+1], numVertices)
  referenceCorners(topologyId, dim, self.baryCenters_[dim+1])

  # compute barycenters
  for codim = 0:dim
    resize!(baryCenters_[codim+1], size(self, codim))
    for i = 1:size(self, codim)
        self.baryCenters_[codim+1][i] = Base.zeros(T,(dim,))
        numCorners = size(self, i, codim, dim)
        for j = 1:numCorners
            baryCenters_[codim+1][i] += baryCenters_[dim+1][subEntity(self, i, codim, j, dim)]
            baryCenters_[codim+1][i] *= 1::T / T( numCorners )
        end
    end
  end

  # compute reference element volume
  self.volume_ = referenceVolume(T, topologyId, dim)

  # compute integration outer normals
  if dim > 0
    resize!(integrationNormals_, size(self, 1))
    referenceIntegrationOuterNormals(topologyId, dim, integrationNormals_)
  end

  # set up geometries
  for i = 1:dim
    createGeometries(self, i, self.geometries_)
  end
end