module DuneGeometry

using Reexport

abstract type AbstractGeometry{T<:Real} end

include("type.impl.jl")
include("type.jl")

include("referenceelement.impl.jl")
include("referenceelement.jl")

# geometries
include("utility.jl")
include("geometry.jl")
include("multilineargeometry.impl.jl")
include("multilineargeometry.jl")
include("affinegeometry.jl")

end # module DuneGeometry
