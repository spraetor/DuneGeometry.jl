using ArgCheck

@enum TopologyConstruction begin
  pyramidConstruction = 0
  prismConstruction = 1
end

"""
  numTopologies(dim)

Obtain the number of topologies of a given dimension `dim`.

**Note:** Valid topology ids are `0,...,numTopologies(dim)-1`.

# Arguments
- `dim::Int`: Dimension.
"""
numTopologies(dim::Int) = (1 << dim)::UInt


"""
  isPyramid(topologyId, dim[, codim=0])

Check whether a pyramid construction was used to create a given codimension `codim`.

# Arguments
- `topologyId::UInt`: id of the topology
- `dim::Int`: Dimension of the topology
- `codim::Int`: Codimension for which the information is desired (defaults to 0)
"""
function isPyramid(topologyId::UInt, dim::Int, codim::Int = 0)
  @argcheck (dim > 0) && (topologyId < numTopologies(dim))
  @argcheck (0 <= codim) && (codim < dim)
 ((topologyId & ~1) & (1u << (dim-codim-1))) == 0
end


"""
  isPrism(topologyId, dim[, codim=0])

Check whether a prism construction was used to create a given codimension `codim`.

# Arguments
- `topologyId::UInt`: id of the topology
- `dim::Int`: Dimension of the topology
- `codim::Int`: Codimension for which the information is desired (defaults to 0)
"""
function isPrism(topologyId::UInt, dim::Int, codim::Int = 0)
  @argcheck (dim > 0) && (topologyId < numTopologies( dim ))
  @argcheck (0 <= codim) && (codim < dim)
  ((topologyId | 1) & (1u << (dim-codim-1))) != 0
end


"""
  baseTopologyId(topologyId, dim[, codim=1])

Obtain the base topology of a given codimension `codim`.

# Arguments
- `topologyId::UInt`: id of the topology
- `dim::Int`: Dimension of the topology
- `codim::Int`: Codimension for which the information is desired (defaults to 0)
"""
function baseTopologyId(topologyId::UInt, dim::Int, codim::Int = 1)
  @argcheck (dim >= 0) && (topologyId < numTopologies( dim ))
  @argcheck (0 <= codim) && (codim <= dim)
  topologyId & ((1u << (dim-codim)) - 1)
end