module DuneGeometry

using Reexport

include("type.impl.jl")
include("type.jl")
@reexport using .Types

include("geometry.jl")
@reexport using .Geometries

include("referenceelement.impl.jl")
include("referenceelement.jl")
@reexport using .ReferenceElements

include("affinegeometry.jl")
@reexport using .AffineGeometries

end # module DuneGeometry
