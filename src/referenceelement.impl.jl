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

@enum TopologyConstruction begin
  pyramidConstruction = 0
  prismConstruction = 1
end


function size(topologyId::UInt, dim::Int, codim::Int)
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


function subTopologyId(topologyId::UInt, dim::Int, codim::Int, i::UInt)
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