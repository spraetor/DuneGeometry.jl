module TypesImpl

using ArgCheck

@enum TopologyConstruction begin
  pyramidConstruction = 0
  prismConstruction = 1
end


# Obtain the number of topologies of a given dimension `dim`.
# Note: Valid topology ids are `0,...,numTopologies(dim)-1`.
numTopologies(dim::Integer) = (1 << dim)


# Check whether a pyramid construction was used to create a given codimension `codim`.
#
# Arguments
# - `topologyId::UInt32`: id of the topology
# - `dim::Integer`: Dimension of the topology
# - `codim::Integer`: Codimension for which the information is desired (defaults to 0)
function isPyramid(topologyId::UInt32, dim::Integer, codim::Integer = 0)::Bool
  @argcheck (dim > 0) && (topologyId < numTopologies(dim))
  @argcheck 0 <= codim < dim
  ((topologyId & ~1) & (1 << (dim-codim-1))) == 0
end


# Check whether a prism construction was used to create a given codimension `codim`.
#
# Arguments
# - `topologyId::UInt32`: id of the topology
# - `dim::Integer`: Dimension of the topology
# - `codim::Integer`: Codimension for which the information is desired (defaults to 0)
function isPrism(topologyId::UInt32, dim::Integer, codim::Integer = 0)::Bool
  @argcheck (dim > 0) && (topologyId < numTopologies(dim))
  @argcheck 0 <= codim < dim
  ((topologyId | 1) & (1 << (dim-codim-1))) != 0
end


# Obtain the base topology of a given codimension `codim`.
#
# Arguments
# - `topologyId::UInt32`: id of the topology
# - `dim::Integer`: Dimension of the topology
# - `codim::Integer`: Codimension for which the information is desired (defaults to 0)
function baseTopologyId(topologyId::UInt32, dim::Integer, codim::Integer = 1)::UInt32
  @argcheck (dim >= 0) && (topologyId < numTopologies(dim))
  @argcheck (0 <= codim) && (codim <= dim)
  topologyId & ((1 << (dim-codim)) - 1)
end

end # module TypesImpl
