struct ReferenceElement{T<:AbstractFloat}
    geometryType::GeometryType
    coordinates::Vector{Vector{T}}
    facets::Vector{Vector{Int}}
    volume::T

    function ReferenceElement(filename::String)
        open(filename) do f
            import YAML
            data = YAML.load_file(f)
            dim = Int(data["dim"])
            new(GeometryType(toBasicType(data["type"]), dim),
                Float64.(data["coordinates"]),
                Float64.(data["facets"])
                Float64(data["volume"]))
        end
    end
end


"Obtain the type of subentity (i,c)"
function type(r::ReferenceElement, i::Int, c::Int)
    dim = dim(r.geometryType)
    if isSimplex(r.geometryType)
        GeometryType(BasicType.simplex, dim - c)
    elseif isCube(r.geometryType)
        GeometryType(BasicType.cube, dim - c)
    else
        if c == 1
            GeometryType(length(r.facets[i]), dim - 1)
        else
            GeometryType(dim - c, dim - c)
        end
    end
end

"Obtain the type of this reference element"
type(r::ReferenceElement) = r.geometryType


"""
    size(g, c)

Number of subentities of codimension `c` in the GeometryType `g`.

# Arguments
- `g::GeometryType`: The type of the geometry.
- `c::Int`: Codimension whose size is requested.
"""
size(g::GeometryType, c::Int) = geometryTypeSizes[basicType(g)][g.dim][c]


"""
    size(r, c)

Number of subentities of codimension `c` in the ReferenceElement `r`.

# Arguments
- `r::ReferenceElement`: The reference element.
- `c::Int`: Codimension whose size is requested.
"""
size(r::ReferenceElement, c::Int) = size(r.geometryType, c)


"""
    size(r, i, c, cc)

Number of subentities of codimension `cc` of subentity `(i,c)`.

Denote by `E` the `i`-th subentity of codimension `c` of the given
reference element `r`. This method returns the number of subentities
of codimension `cc` of the reference element, that are also
a subentity of `E`. If `cc<c` this number is zero.

# Arguments
- `r::ReferenceElement`: The reference element.
- `i:Int`: The number of subentity `E` (`0 <= i < size(r,c)`)
- `c::Int`: Codimension of subentity `E` (`0 <= c <= dim(r)`)
- `cc::Int`: Codimension whose size is desired (`0 <= cc <= dim(r)`)
"""
size(r::ReferenceElement, i::Int, c::Int, cc::Int) = size(type(r, i, c), cc - c)


"""
    subEntity(i,c,ii,cc)

Obtain number of `ii`-th subentity with codim `cc` of `(i,c)`.

Denote by E the i-th subentity of codimension c of the current
reference element. And denote by S the ii-th subentity of codimension
(cc-c) of E. Then, S is a also a subentity of codimension cc of the current
reference element. This method returns the number of S with respect
to the current reference element.

# Arguments
- `r::ReferenceElement`: The reference element.
- `i::Int`: Number of subentity E (0 <= i < size( c ))
- `c::Int`: Codimension of subentity E
- `ii:Int`: Number of subentity S (with respect to E)
- `cc:Int`: Codimension of subentity S (c <= cc <= dim)
"""
subEntity(r::ReferenceElement, i::Int, c::Int, ii::Int, cc::Int)




function initialize{T<:Real}(r::ReferenceElement{T}, topologyId::UInt)
  @argcheck topologyId < numTopologies(r.dim)

  # set up subentities
  for codim = 0:r.dim
    s = size(topologyId, dim, codim)
    resize!(r.info_[codim+1], s)

    for i = 1:s
      info_[codim+1][i].initialize(topologyId, codim, i)
    end
  end

  # compute corners
  numVertices = size(r,dim)
  resize!(r.baryCenters_[r.dim+1], numVertices)
  referenceCorners(topologyId, dim, r.baryCenters_[r.dim+1])

  # compute barycenters
  for codim = 0:r.dim
    resize!(baryCenters_[codim+1], size(r, codim))
    for i = 1:size(r, codim)
        r.baryCenters_[codim+1][i] = Base.zeros(T,(r.dim,))
        numCorners = size(r, i, codim, dim)
        for j = 1:numCorners
            baryCenters_[codim+1][i] += baryCenters_[r.dim+1][subEntity(i, codim, j, dim)]
            baryCenters_[codim+1][i] *= 1::T / T( numCorners )
        end
    end
  end

  # compute reference element volume
  r.volume_ = referenceVolume(T, topologyId, dim)

  # compute integration outer normals
  if dim > 0
    resize!(integrationNormals_, size(r, 1))
    referenceIntegrationOuterNormals(topologyId, dim, integrationNormals_)
  end

  # set up geometries
  for i = 1:r.dim
    createGeometries(r, i, r.geometries_)
  end
end