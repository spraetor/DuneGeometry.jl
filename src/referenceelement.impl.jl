using ArgCheck

# internal: Number of sub-entities of each co-dimension
# [BasicType][dim][codim]
geometryTypeSizes = Dict(
    BasicType.simplex => [
        [1],            # dim == 0
        [1, 2],         # dim == 1
        [1, 3, 3],      # dim == 2
        [1, 4, 6, 4]],  # dim == 3
    BasicType.cube => [
        [1],            # dim == 0
        [1, 2],         # dim == 1
        [1, 4, 4],      # dim == 2
        [1, 6, 12, 8]], # dim == 3
    BasicType.prism => [
        [], [], [],
        [1, 5, 9, 6]],  # dim == 3
    BasicType.pyramid => [
        [], [], [],
        [1, 5, 8, 5]]   # dim == 3
)


function size(topologyId::UInt, dim::Int, codim::Int)::UInt
  @argcheck (dim >= 0) && (topologyId < numTopologies(dim))
  @argcheck (0 <= codim) && (codim <= dim)
  
  if codim > 0
    baseId = baseTopologyId(topologyId, dim)
    m = size(baseId, dim-1, codim-1)
    if isPrism(topologyId, dim)
      n = (codim < dim ? size( baseId, dim-1, codim ) : 0)
      (n + 2*m)::UInt
    else
      @assert isPyramid(topologyId, dim)
      n = (codim < dim ? size( baseId, dim-1, codim ) : 1)
      (m + n)::UInt
    end
  else
    return 1::UInt
  end
end


function subTopologyId(topologyId::UInt, dim::Int, codim::Int, i::UInt)::UInt
  @argcheck i < size(topologyId, dim, codim)
  mydim = dim - codim

  if codim > 0
    baseId = baseTopologyId(topologyId, dim)
    m = size(baseId, dim-1, codim-1)

    if isPrism(topologyId, dim)
      n = (codim < dim ? size(baseId, dim-1, codim) : 0)
      if i < n
        return subTopologyId(baseId, dim-1, codim, i) | (prismConstruction::UInt << (mydim - 1))
      else
        return subTopologyId(baseId, dim-1, codim-1, (i < n+m ? i-n : i-(n+m)) )
      end
    else
      @assert isPyramid(topologyId, dim)
      if i < m
        return subTopologyId(baseId, dim-1, codim-1, i)
      elseif codim < dim
        return subTopologyId(baseId, dim-1, codim, i-m) | (pyramidConstruction::UInt << (mydim - 1))
      else
        return 0::UInt;
      end
    end
  else
    return topologyId
  end
end


function subTopologyNumbering(topologyId::UInt, dim::Int, codim::Int, i::UInt, subcodim::Int, out::AbstractVector{UInt})
  @argcheck (codim >= 0) && (subcodim >= 0) && (codim + subcodim <= dim)
  @argcheck i < size(topologyId, dim, codim)
  @argcheck length(out) == size(subTopologyId(topologyId, dim, codim, i), dim-codim, subcodim)

  if codim == 0
    out = range(0,length=length(out))
  elseif subcodim == 0
    @assert length(out) == 1
    out[1] = i
  else
    beginOut = 1
    endOut = length(out)
    baseId = baseTopologyId(topologyId, dim)
    m = size(baseId, dim-1, codim-1)
    mb = size(baseId, dim-1, codim+subcodim-1)
    nb = (codim + subcodim < dim ? size(baseId, dim-1, codim+subcodim) : 0)

    if isPrism(topologyId, dim)
      n = size(baseId, dim-1, codim)
      if i < n
        subId = subTopologyId(baseId, dim-1, codim, i)

        beginBase = beginOut
        if codim + subcodim < dim
          beginBase = beginOut + size(subId, dim-codim-1, subcodim)
          subTopologyNumbering(baseId, dim-1, codim, i, subcodim, view(out,beginOut:beginBase))
        end

        ms = size(subId, dim-codim-1, subcodim-1)
        subTopologyNumbering(baseId, dim-1, codim, i, subcodim-1, view(out,beginBase:beginBase+ms))
        for j = 1:ms
          out[beginBase + j - 1] += nb
          out[beginBase + j + ms - 1] = out[beginBase + j -1] + mb
        end
      else  
        s = (i < n+m ? 0 : 1)
        subTopologyNumbering(baseId, dim-1, codim-1, i-(n+s*m), subcodim, out)
        for j in eachIndex(out)
          out[j] += nb + s*mb
        end
      end
    else
      @assert isPyramid(topologyId, dim)

      if i < m
        subTopologyNumbering(baseId, dim-1, codim-1, i, subcodim, out)
      else
        ubId = subTopologyId(baseId, dim-1, codim, i-m)
        ms = size(subId, dim-codim-1, subcodim-1)

        subTopologyNumbering(baseId, dim-1, codim, i-m, subcodim-1, view(out,beginOut:beginOut+ms))
        if codim+subcodim < dim
          subTopologyNumbering(baseId, dim-1, codim, i-m, subcodim, view(out,beginOut+ms:endOut))
          for j = (beginOut + ms):endOut
            out[j] += mb
          end      
        else
          out[beginOut + ms] = mb
        end
      end
    end
  end
end


function checkInside{T<:Real}(topologyId::UInt, dim::Int, x::AbstractVector{T}, tolerance::T, factor::T = 1::T)::Bool
  @argcheck (dim >= 0) && (dim <= cdim)
  @argcheck topologyId < numTopologies(dim)

  if dim > 0
    baseFactor = (isPrism(topologyId, dim) ? factor : factor - x[dim])
    if (x[dim] > -tolerance) && (factor - x[dim] > -tolerance)
      return checkInside(baseTopologyId(topologyId, dim), dim-1, x, tolerance, baseFactor)
    else
      return false
    end
  else
    return true
  end
end


function referenceCorners{T<:Real}(topologyId::UInt, dim::Int, corners::AbstractVector{AbstractArray{T,1}})::UInt
  @argcheck (dim >= 0) && (dim <= length(corners[1]))
  @argcheck topologyId < numTopologies(dim)

  if dim > 0
    nBaseCorners = referenceCorners(baseTopologyId(topologyId, dim), dim-1, corners)
    @assert nBaseCorners == size(baseTopologyId(topologyId, dim), dim-1, dim-1)

    if isPrism(topologyId, dim)
      for j = 1:nBaseCorners
        corners[nBaseCorners+j] = corners[j]
        corners[nBaseCorners+j][dim] = 1::T
      end
      return 2*nBaseCorners
    else
      corners[nBaseCorners+1] = zeros(T,Base.size(corners[nBaseCorners+1]))
      corners[nBaseCorners+1][dim] = 1::T
      return nBaseCorners+1
    end
  else
    corners[1] = zeros(T,Base.size(corners[1]))
    return 1
  end
end


function referenceVolumeInverse(topologyId::UInt, dim::Int)::UInt
  @argcheck (dim >= 0) && (topologyId < numTopologies(dim))

  if dim > 0
      baseValue = referenceVolumeInverse(baseTopologyId(topologyId, dim), dim-1)
      return (isPrism(topologyId, dim) ? baseValue : baseValue * UInt(dim))
  else
    return 1
  end
end

referenceVolume(::Type{T}, topologyId::UInt, dim::Int)::T where {T<:Real} = 1::T / T(referenceVolumeInverse(topologyId,dim))


function referenceOrigins{T<:Real}(topologyId::UInt, dim::Int, codim::Int, origins::AbstractVector{AbstractArray{T,1}})::UInt
  cdim = length(origins[1])
  @argcheck (dim >= 0) && (dim <= cdim)
  @argcheck topologyId < numTopologies( dim )
  @argcheck (codim >= 0) && (codim <= dim)

  if codim > 0
    baseId = baseTopologyId(topologyId, dim)
    if isPrism(topologyId, dim)
      n = (codim < dim ? referenceOrigins(baseId, dim-1, codim, origins) : 0)
      m = referenceOrigins(baseId, dim-1, codim-1, view(origins,n+1:length(origins)))
      for j = 1:m
        origins[n+m+j] = origins[n+j]
        origins[n+m+j][dim] = 1::T
      end
      return n+2*m
    else
      m = referenceOrigins(baseId, dim-1, codim-1, origins)
      if codim == dim
        origins[m+1] = zeros(T,Base.size(origins[m+1]))
        origins[m+1][dim] = 1::T
        return m+1
      else
        return m + referenceOrigins(baseId, dim-1, codim, view(origins,m+1:length(origins)))
      end
    end
  else
    origins[1] = zeros(T,Base.size(origins[1]))
    return 1
  end
end


function referenceEmbeddings{t<:Real}(topologyId::UInt, dim::Int, codim::Int, 
                                      origins::AbstractVector{AbstractArray{T,1}},
                                      jacobianTransposeds::AbstractVector{AbstractArray{T,2}})::UInt
  cdim = length(origins[1])
  mydim = Base.size(jacobianTransposeds,1)
  @argcheck cdim == Base.size(jacobianTransposeds,2)
  @argcheck (0 <= codim) && (codim <= dim) && (dim <= cdim)
  @argcheck (dim - codim <= mydim) && (mydim <= cdim)
  @argcheck topologyId < numTopologies(dim)

  if (0 < codim) && (codim <= dim)
    baseId = baseTopologyId(topologyId, dim)
    if isPrism( topologyId, dim )
      n = (codim < dim ? referenceEmbeddings(baseId, dim-1, codim, origins, jacobianTransposeds) : 0)
      jacobianTransposeds[1:n][dim-codim,dim] = 1::T

      m = referenceEmbeddings(baseId, dim-1, codim-1, 
        view(origins,(n+1):length(origins)), 
        view(jacobianTransposeds,(n+1):length(jacobianTransposeds)))

      for j = 1:m
        origins[n+m+j] = origins[n+j]
        origins[n+m+j][dim] = 1::T
        jacobianTransposeds[n+m+j] = jacobianTransposeds[n+j]
      end
      return n+2*m
    else # !isPrism
      m = referenceEmbeddings(baseId, dim-1, codim-1, origins, jacobianTransposeds)
      if codim == dim
        origins[m+1] = zeros(T, Base.size(origins[m+1]))
        origins[m+1][dim] = 1::T
        jacobianTransposeds[m+1] = zeros(T, Base.size(jacobianTransposeds[m+1]))
        return m+1
      elseif codim < dim
        n = referenceEmbeddings(baseId, dim-1, codim, 
          view(origins,(m+1):length(origins)), 
          view(jacobianTransposeds,(m+1):length(jacobianTransposeds)))
        for j = 1:n
          jacobianTransposeds[m+j][dim-codim,:] = -origins[m+j]
          jacobianTransposeds[m+j][dim-codim,dim] = 1::T
        end
        return m+n
      end
    end
  elseif codim == 0
    origins[1] = zeros(T,Base.size(origins[1]))
    jacobianTransposeds[1] = zeros(T,Base.size(jacobianTransposeds[1]))
    for k = 1:dim
      jacobianTransposeds[1][k,k] = 1::T
    end
    return 1
  end

  # this point should not be reached since all cases are handled before.
  return 0
end


function referenceIntegrationOuterNormals{T<:Real}(topologyId::UInt, dim::Int, 
                                                   origins::AbstractVector{AbstractArray{T,1}},
                                                   normals::AbstractVector{AbstractArray{T,1}})::UInt
  cdim = length(origins[1])
  @argcheck cdim == length(normals[1])
  @argcheck (dim > 0) && (dim <= cdim)
  @argcheck topologyId < numTopologies(dim)

  if dim > 1
    baseId = baseTopologyId(topologyId, dim)
    if isPrism(topologyId, dim)
      numBaseFaces = referenceIntegrationOuterNormals(baseId, dim-1, origins, normals)
      for i = 1:2
        normals[numBaseFaces+i] = zeros(T,Base.size(normals[numBaseFaces+i+1]))
        normals[numBaseFaces+i][dim] = T(2*(i-1) - 1)
      end
      return numBaseFaces+2
    else
      normals[1] = zeros(T,Base.size(normals[1]))
      normals[1][dim] = -1::T;
      numBaseFaces = referenceIntegrationOuterNormals(baseId, dim-1, 
        view(origins,2:length(origins)), view(normals,2:length(normals)))
      for i = 1:numBaseFaces
        normals[i+1][dim] = dot(normals[i+1], origins[i+1])
      end
      return numBaseFaces+1
    end
  else
    for i = 1:2
      normals[i] = zeros(T,Base.size(normals[i+1]))
      normals[i][1] = T(2*(i-1) - 1)
    end
    return 2
  end
end


function referenceIntegrationOuterNormals{T<:Real}(topologyId::UInt, dim::Int,
                                                   normals::AbstractVector{AbstractArray{T,1}})::UInt
  cdim = length(normals[1])                                     
  @argcheck (dim > 0) && (dim <= cdim)

  origins = [Vector{T}(undef, cdim) for _ in 1:size(topologyId, dim, 1)]
  referenceOrigins(topologyId, dim, 1, origins)

  numFaces = referenceIntegrationOuterNormals(topologyId, dim, origins, normals)
  @assert numFaces == size(topologyId, dim, 1)

  return numFaces
end