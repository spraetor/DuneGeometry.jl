module ReferenceElementsImpl

using ..Types: GeometryType
import ..TypesImpl

using ArgCheck

function size(topologyId::UInt32, dim::Integer, codim::Integer)
  @argcheck (dim >= 0) && (topologyId < TypesImpl.numTopologies(dim))
  @argcheck (0 <= codim) && (codim <= dim)

  if codim > 0
    baseId = TypesImpl.baseTopologyId(topologyId, dim)
    m = size(baseId, dim-1, codim-1)
    if TypesImpl.isPrism(topologyId, dim)
      n = (codim < dim ? size( baseId, dim-1, codim ) : 0)
      return n + 2*m
    else
      @assert TypesImpl.isPyramid(topologyId, dim)
      n = (codim < dim ? size( baseId, dim-1, codim ) : 1)
      return m + n
    end
  else
    return 1
  end
end


function subTopologyId(topologyId::UInt32, dim::Integer, codim::Integer, i::Integer)
  @argcheck i < size(topologyId, dim, codim)
  mydim = dim - codim

  if codim > 0
    baseId = TypesImpl.baseTopologyId(topologyId, dim)
    m = size(baseId, dim-1, codim-1)

    if isPrism(topologyId, dim)
      n = (codim < dim ? size(baseId, dim-1, codim) : 0)
      if i < n
        return subTopologyId(baseId, dim-1, codim, i) | (prismConstruction::UInt << (mydim - 1))
      else
        return subTopologyId(baseId, dim-1, codim-1, (i < n+m ? i-n : i-(n+m)) )
      end
    else
      @assert TypesImpl.isPyramid(topologyId, dim)
      if i < m
        return subTopologyId(baseId, dim-1, codim-1, i)
      elseif codim < dim
        return subTopologyId(baseId, dim-1, codim, i-m) | (pyramidConstruction::UInt << (mydim - 1))
      else
        return 0;
      end
    end
  else
    return topologyId
  end
end


function subTopologyNumbering!(topologyId::UInt32, dim::Integer, codim::Integer, i::Integer,
                               subcodim::Integer,
                   #= inout =# out::AbstractVector{UInt})
  @argcheck (codim >= 0) && (subcodim >= 0) && (codim + subcodim <= dim)
  @argcheck i <= size(topologyId, dim, codim)
  @argcheck length(out) == size(subTopologyId(topologyId, dim, codim, i), dim-codim, subcodim)

  if codim == 0
    out = Base.range(0,length=length(out))
  elseif subcodim == 0
    @assert length(out) == 1
    out[1] = i
  else
    beginOut = 1
    endOut = length(out)
    baseId = TypesImpl.baseTopologyId(topologyId, dim)
    m = size(baseId, dim-1, codim-1)
    mb = size(baseId, dim-1, codim+subcodim-1)
    nb = (codim + subcodim < dim ? size(baseId, dim-1, codim+subcodim) : 0)

    if TypesImpl.isPrism(topologyId, dim)
      n = size(baseId, dim-1, codim)
      if i < n
        subId = subTopologyId(baseId, dim-1, codim, i)

        beginBase = beginOut
        if codim + subcodim < dim
          beginBase = beginOut + size(subId, dim-codim-1, subcodim)
          subTopologyNumbering!(baseId, dim-1, codim, i, subcodim,
            Base.view(out,beginOut:beginBase))
        end

        ms = size(subId, dim-codim-1, subcodim-1)
        subTopologyNumbering!(baseId, dim-1, codim, i, subcodim-1,
          Base.view(out,beginBase:beginBase+ms))
        for j = 1:ms
          out[beginBase + j - 1] += nb
          out[beginBase + j + ms - 1] = out[beginBase + j -1] + mb
        end
      else
        s = (i < n+m ? 0 : 1)
        subTopologyNumbering!(baseId, dim-1, codim-1, i-(n+s*m), subcodim, out)
        for j in eachIndex(out)
          out[j] += nb + s*mb
        end
      end
    else
      @assert TypesImpl.isPyramid(topologyId, dim)

      if i < m
        subTopologyNumbering!(baseId, dim-1, codim-1, i, subcodim, out)
      else
        subId = subTopologyId(baseId, dim-1, codim, i-m)
        ms = size(subId, dim-codim-1, subcodim-1)

        subTopologyNumbering!(baseId, dim-1, codim, i-m, subcodim-1,
          Base.view(out,beginOut:beginOut+ms))
        if codim+subcodim < dim
          subTopologyNumbering!(baseId, dim-1, codim, i-m, subcodim,
            Base.view(out,beginOut+ms:endOut))
          for j = (beginOut + ms):endOut
            out[j] += mb
          end
        else
          out[beginOut + ms] = mb
        end
      end
    end
  end
  return nothing
end


function checkInside(topologyId::UInt32, dim::Integer, x::AbstractVector{Real}, tolerance::Real, factor::Real = 1)::Bool
  @argcheck 0 <= dim <= cdim
  @argcheck topologyId < TypesImpl.numTopologies(dim)

  if dim > 0
    baseFactor = (TypesImpl.isPrism(topologyId, dim) ? factor : factor - x[dim])
    if (x[dim] > -tolerance) && (factor - x[dim] > -tolerance)
      return checkInside(TypesImpl.baseTopologyId(topologyId, dim), dim-1, x, tolerance, baseFactor)
    else
      return false
    end
  else
    return true
  end
end


function referenceCorners!(topologyId::UInt32, dim::Integer,
               #= inout =# corners::AbstractVector{AbstractArray{Real,1}})
  @argcheck 0 <= dim <= length(corners[1])
  @argcheck topologyId < TypesImpl.numTopologies(dim)

  if dim > 0
    nBaseCorners = referenceCorners!(TypesImpl.baseTopologyId(topologyId, dim), dim-1, corners)
    @assert nBaseCorners == size(TypesImpl.baseTopologyId(topologyId, dim), dim-1, dim-1)

    if TypesImpl.isPrism(topologyId, dim)
      for j = 1:nBaseCorners
        corners[nBaseCorners+j] = copy(corners[j])
        corners[nBaseCorners+j][dim] = 1
      end
      return 2*nBaseCorners
    else
      corners[nBaseCorners+1] .= 0
      corners[nBaseCorners+1][dim] = 1
      return nBaseCorners+1
    end
  else
    corners[1] .= 0
    return 1
  end
end


function referenceVolumeInverse(topologyId::UInt32, dim::Integer)
  @argcheck (dim >= 0) && (topologyId < TypesImpl.numTopologies(dim))

  if dim > 0
      baseValue = referenceVolumeInverse(TypesImpl.baseTopologyId(topologyId, dim), dim-1)
      return (TypesImpl.isPrism(topologyId, dim) ? baseValue : baseValue * dim)
  else
    return 1
  end
end

function referenceVolume(::Type{T}, topologyId::UInt32, dim::Integer)::T where {T<:Real}
  1::T / T(referenceVolumeInverse(topologyId,dim))
end


function referenceOrigins!(topologyId::UInt32, dim::Integer, codim::Integer,
                #= inout=# origins::AbstractVector{AbstractArray{Real,1}})
  cdim = length(origins[1])
  @argcheck 0 <= dim <= cdim
  @argcheck topologyId < TypesImpl.numTopologies( dim )
  @argcheck 0 <= codim <= dim

  if codim > 0
    baseId = TypesImpl.baseTopologyId(topologyId, dim)
    if TypesImpl.isPrism(topologyId, dim)
      n = (codim < dim ? referenceOrigins!(baseId, dim-1, codim, origins) : 0)
      m = referenceOrigins!(baseId, dim-1, codim-1, Base.view(origins,n+1:length(origins)))
      for j = 1:m
        origins[n+m+j] = copy(origins[n+j])
        origins[n+m+j][dim] = 1
      end
      return n+2*m
    else
      m = referenceOrigins!(baseId, dim-1, codim-1, origins)
      if codim == dim
        origins[m+1] .= 0
        origins[m+1][dim] = 1
        return m+1
      else
        return m + referenceOrigins!(baseId, dim-1, codim, Base.view(origins,m+1:length(origins)))
      end
    end
  else
    origins[1] .= 0
    return 1
  end
end


function referenceEmbeddings!(topologyId::UInt32, dim::Integer, codim::Integer,
                              origins::AbstractVector{AbstractArray{Real,1}},
                              jacobianTransposeds::AbstractVector{AbstractArray{Real,2}})
  cdim = length(origins[1])
  mydim = Base.size(jacobianTransposeds,1)
  @argcheck cdim == Base.size(jacobianTransposeds,2)
  @argcheck 0 <= codim <= dim <= cdim
  @argcheck dim - codim <= mydim <= cdim
  @argcheck topologyId < TypesImpl.numTopologies(dim)

  if 0 < codim <= dim
    baseId = TypesImpl.baseTopologyId(topologyId, dim)
    if TypesImpl.isPrism( topologyId, dim )
      n = (codim < dim ? referenceEmbeddings!(baseId, dim-1, codim, origins, jacobianTransposeds) : 0)
      jacobianTransposeds[1:n][dim-codim,dim] .= 1

      m = referenceEmbeddings!(baseId, dim-1, codim-1,
        Base.view(origins,(n+1):length(origins)),
        Base.view(jacobianTransposeds,(n+1):length(jacobianTransposeds)))

      for j = 1:m
        origins[n+m+j] = copy(origins[n+j])
        origins[n+m+j][dim] = 1
        jacobianTransposeds[n+m+j] = copy(jacobianTransposeds[n+j])
      end
      return n+2*m
    else # !isPrism
      m = referenceEmbeddings!(baseId, dim-1, codim-1, origins, jacobianTransposeds)
      if codim == dim
        origins[m+1] .= 0
        origins[m+1][dim] = 1
        jacobianTransposeds[m+1] .= 0
        return m+1
      elseif codim < dim
        n = referenceEmbeddings!(baseId, dim-1, codim,
          Base.view(origins,(m+1):length(origins)),
          Base.view(jacobianTransposeds,(m+1):length(jacobianTransposeds)))
        for j = 1:n
          jacobianTransposeds[m+j][dim-codim,:] = -origins[m+j]
          jacobianTransposeds[m+j][dim-codim,dim] = 1
        end
        return m+n
      end
    end
  elseif codim == 0
    origins[1] .= 0
    jacobianTransposeds[1] .= 0
    for k = 1:dim
      jacobianTransposeds[1][k,k] = 1
    end
    return 1
  end

  # this point should not be reached since all cases are handled before.
  return 0
end


function referenceIntegrationOuterNormals!(topologyId::UInt32, dim::Integer,
                                  #= in =# origins::AbstractVector{AbstractArray{Real,1}},
                                #= inout=# normals::AbstractVector{AbstractArray{Real,1}})
  cdim = length(origins[1])
  @argcheck cdim == length(normals[1])
  @argcheck 0 < dim <= cdim
  @argcheck topologyId < TypesImpl.numTopologies(dim)

  if dim > 1
    baseId = TypesImpl.baseTopologyId(topologyId, dim)
    if TypesImpl.isPrism(topologyId, dim)
      numBaseFaces = referenceIntegrationOuterNormals!(baseId, dim-1, origins, normals)
      for i = 1:2
        normals[numBaseFaces+i] .= 0
        normals[numBaseFaces+i][dim] = 2*(i-1) - 1
      end
      return numBaseFaces+2
    else
      normals[1] .= 0
      normals[1][dim] = -1
      numBaseFaces = referenceIntegrationOuterNormals!(baseId, dim-1,
        Base.view(origins,2:length(origins)), Base.view(normals,2:length(normals)))
      for i = 1:numBaseFaces
        normals[i+1][dim] = dot(normals[i+1], origins[i+1])
      end
      return numBaseFaces+1
    end
  else
    for i = 1:2
      normals[i] .= 0
      normals[i][1] = 2*(i-1) - 1
    end
    return 2
  end
end


function referenceIntegrationOuterNormals!(topologyId::UInt32, dim::Integer,
                          #= inout=# normals::AbstractVector{AbstractArray{T,1}}) where {T<:Real}
  cdim = length(normals[1])
  @argcheck 0 < dim <= cdim

  origins = [Vector{T}(undef, cdim) for _ in 1:size(topologyId, dim, 1)]
  referenceOrigins!(topologyId, dim, 1, origins)

  numFaces = referenceIntegrationOuterNormals!(topologyId, dim, origins, normals)
  @assert numFaces == size(topologyId, dim, 1)

  return numFaces
end


maxSubEntityCount(dim::Integer) = max(binomial.(dim,0:dim) .* (1 .<< 0:dim))

SubEntityFlags = BitArray{1}

struct SubEntityInfo
  # public
  dimension::Int

  # private
  numbering_::Vector{UInt}
  offset_::Vector{UInt}
  type_::GeometryType
  containsSubentity_::Vector{SubEntityFlags}
  maxSubEntityCount_::Int

  "Constructor"
  SubEntityInfo(dim::Integer) = new(dim,      # dimension
    undef,                                    # numbering_
    Vector{UInt}(undef, dim+2),               # offset_
    undef,                                    # type_
    Vector{SubEntityFlags}(undef, dim+1),     # containsSubentity_
    maxSubEntityCount(dim))                   # maxSubEntityCount_
end


function size(self::SubEntityInfo, cc::Integer)
  @argcheck 0 <= cc <= dim
  return self.offset_[cc+2] - self.offset_[cc+1]
end


function number(self::SubEntityInfo, ii::Integer, cc::Integer)
  @argcheck 0 < ii <= size(self, cc)
  return self.numbering_[self.offset_[cc+1] + ii]
end


function numbers(self::SubEntityInfo, cc::Integer)
  return Base.view(self.numbering_, self.offset_[cc+1]:self.offset_[cc+2]), self.containsSubentity_[cc+1]
end


function initialize!(self::SubEntityInfo, topologyId::UInt32, codim::Integer, i::Integer)
  dim = self.dimension
  subId = subTopologyId(topologyId, dim, codim, i)
  self.type_ = GeometryType(subId, dim-codim)

  # compute offsets
  self.offset_ .= 0
  for cc = codim:dim
    self.offset_[cc+2] = self.offset_[cc+1] + size(subId, dim-codim, cc-codim)
  end

  # compute subnumbering
  resize!(self.numbering_, self.offset_[dim+2])
  for cc = codim:dim
    subTopologyNumbering(topologyId, dim, codim, i, cc-codim,
      Base.view(self.numbering_, self.offset_[cc+1]:self.offset_[cc+2]))
  end

  # initialize containsSubentity lookup-table
  for cc = 0:dim
    self.containsSubentity_[cc+1] = Base.falses(self.maxSubEntityCount_)
    self.containsSubentity_[cc+1][number.(self,1:size(self,cc),cc)] .= true
  end
  return nothing
end

end # module ReferenceElementsImpl
