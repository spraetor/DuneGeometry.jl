module DuneGeometry

using Reexport

abstract type AbstractGeometry{T<:Real} end

include("type.impl.jl")
include("type.jl")

include("referenceelement.impl.jl")
include("referenceelement.jl")

# geometries
include("affinegeometry.jl")
include("geometry.jl")

end # module DuneGeometry
