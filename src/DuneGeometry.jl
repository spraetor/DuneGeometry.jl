module DuneGeometry

export ReferenceElement,GeometryType,BasicType

include("type.impl.jl")
include("type.jl")
using .Types

include("geometry.jl")
using .Geometries

include("referenceelement.impl.jl")
include("referenceelement.jl")
using .ReferenceElements

include("affinegeometry.jl")
using .AffineGeometries

end # module DuneGeometry
