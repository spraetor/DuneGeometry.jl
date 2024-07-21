# types
export ReferenceElement
# methods
export size,subEntity,subEntities,type,position,checkInside,geometry,volume,integrationOuterNormal

using ArgCheck

# Import some implementation details
import ..ReferenceElementsImpl
import ..TypesImpl

# Some type aliases
Coordinate{T<:Real} = Vector{T}


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

    # private members
    volume_::T
    baryCenters_::Vector{Vector{Coordinate{T}}}
    integrationNormals_::Vector{Coordinate{T}}
    info_::Vector{Vector{ReferenceElementsImpl.SubEntityInfo}}

    # geometry information
    origins_::Vector{Vector{Coordinate{T}}}
    jacobianTransposeds_::Vector{Vector{Array{T,2}}}
end


"""
    size(r, c)

Number of subentities of codimension `c` in the ReferenceElement `r`.

# Arguments
- `r::ReferenceElement`: The reference element.
- `c::Int`: Codimension whose size is requested.
"""
function size(self::ReferenceElement{T}, c::Integer) where {T<:Real}
    @argcheck 0 <= c <= self.dimension
    return length(self.info_[c+1])
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
function size(self::ReferenceElement{T}, i::Integer, c::Integer, cc::Integer) where {T<:Real}
    @argcheck 0 < i <= size(self, c)
    return ReferenceElementsImpl.size(self.info_[c+1][i], cc)
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
function subEntity(self::ReferenceElement{T}, i::Integer, c::Integer, ii::Integer, cc::Integer) where {T<:Real}
    @argcheck 0 < i <= size(self, c)
    return ReferenceElementsImpl.number(self.info_[c+1][i], ii, cc )
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
function subEntities(self::ReferenceElement{T}, i::Integer, c::Integer, cc::Integer) where {T<:Real}
    @argcheck 0 < i <= size(self, c)
    return ReferenceElementsImpl.numbers(self.info_[c+1][i], cc)
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
function type(self::ReferenceElement{T}, i::Integer, c::Integer)::GeometryType where {T<:Real}
    @argcheck 0 < i <= size(self, c)
    return self.info_[c+1][i].type
end


"Obtain the type of this reference element"
function type(self::ReferenceElement{T})::GeometryType where {T<:Real}
    return type(self, 1, 0)
end


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
function position(self::ReferenceElement{T}, i::Integer, c::Integer)::Vector{T} where {T<:Real}
    @argcheck 0 <= c <= self.dimension
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
function checkInside(self::ReferenceElement{T}, x::DenseVector{T})::Bool where {T<:Real}
    return ReferenceElementsImpl.checkInside(type(self).topologyId, self.dimension, x, 64 * eps(T))
end


"""
    geometry{G}(r, i)

Obtain the embedding of subentity `(i,codim)` into the reference element

Denote by E the i-th subentity of codimension codim of the current
reference element. This method returns a geometry of type `G`
that maps the reference element of E into the current reference element.

# Arguments
- `i::Int`: number of subentity E (0 < i <= size( c ))
- `c::Int`: Codimension of subentity E
"""
function geometry(::Type{G}, self::ReferenceElement{T}, i::Integer, c::Integer)::G where {T<:Real,G<:AbstractGeometry{T}}
    return G(subRefElement(self,i,c), self.origins_[c+1][i], self.jacobianTransposeds_[c+1][i])
end


"Obtain the volume of the reference element."
function volume(self::ReferenceElement{T})::T where {T<:Real}
    self.volume_
end


"""
    integrationOuterNormal(r, face)

Obtain the integration outer normal of the reference element.

The integration outer normal is the outer normal whose length coincides
with the face's integration element.

# Arguments
- `face::Int`: Index of the face, whose normal is desired.
"""
function integrationOuterNormal(self::ReferenceElement{T}, face::Integer)::Vector{T} where {T<:Real}
    @argcheck 0 < face <= Base.length(self.integrationNormals_)
    return self.integrationNormals_[face]
end


# constructor of ReferenceElement
function ReferenceElement{T}(type::GeometryType) where {T<:Real}
    dimension = type.dim
    volume = ReferenceElementsImpl.referenceVolume(T, type.topologyId, dimension)

    info = [ Vector{ReferenceElementsImpl.SubEntityInfo}() for _ in 1:dimension+1 ]
    for codim = 0:dimension
        s = ReferenceElementsImpl.size(type.topologyId, dimension, codim)
        resize!(info[codim+1], s)
        for i = 1:s
            info[codim+1][i] = ReferenceElementsImpl.SubEntityInfo(type, codim, i)
        end
    end

    baryCenters = Vector{Vector{Coordinate{T}}}(undef, dimension+1)

    # compute corners
    numVertices = length(info[dimension+1])
    baryCenters[dimension+1] = [Vector{T}(undef,dimension) for _ = 1:numVertices]
    ReferenceElementsImpl.referenceCorners!(type.topologyId, dimension, baryCenters[dimension+1])

    # compute barycenters
    for codim = 0:dimension-1
        baryCenters[codim+1] = [Vector{T}(undef,dimension) for _ = 1:length(info[codim+1])]
        for i = 1:length(info[codim+1])
            baryCenters[codim+1][i] .= 0
            numCorners = ReferenceElementsImpl.size(info[codim+1][i], dimension)
            for j = 1:numCorners
                k = ReferenceElementsImpl.number(info[codim+1][i], j, dimension)
                baryCenters[codim+1][i] += baryCenters[dimension+1][k]
            end
            baryCenters[codim+1][i] *= T(1) / T( numCorners )
        end
    end

    # compute integration outer normals
    integrationNormals = Vector{Coordinate{T}}()
    if dimension > 0
      resize!(integrationNormals, length(info[2]))
      fill!(integrationNormals, Vector{T}(undef,dimension))
      ReferenceElementsImpl.referenceIntegrationOuterNormals!(type.topologyId, dimension, integrationNormals)
    end

    origins = Vector{Vector{Coordinate{T}}}(undef, dimension+1)
    jacobianTransposeds = Vector{Vector{Array{T,2}}}(undef, dimension+1)

    ref = ReferenceElement{T}(dimension, volume, baryCenters, integrationNormals, info, origins, jacobianTransposeds)

    for codim = 0:dimension
        createGeometries!(ref,codim)
    end

    return ref
end


function subRefElement(self::ReferenceElement{T}, i::Integer, cc::Integer) where {T<:Real}
    return cc == 0 ? self : ReferenceElement{T}(type(self,i,cc))
end

function createGeometries!(self::ReferenceElement{T}, codim::Integer) where {T<:Real}
    dim = self.dimension
    s = size(self, codim)
    self.origins_[codim+1] = [ Vector{T}(undef,dim) for _ = 1:s ]
    self.jacobianTransposeds_[codim+1] = [ Matrix{T}(undef, dim-codim, dim) for _ in 1:s ]
    ReferenceElementsImpl.referenceEmbeddings!(type(self).topologyId, dim, codim, self.origins_[codim+1], self.jacobianTransposeds_[codim+1])
end
