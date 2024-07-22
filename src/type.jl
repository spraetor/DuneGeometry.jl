# types
export BasicType,GeometryType
# methods
export basicType,toId,isVertex,isLine,isTriangle,isQuadrilateral,isTetrahedron,isHexahedron,isPyramid,isPrism,isSimplex,isCube,isConical,isPrismatic,isNone

using EnumX

# Each entity can be tagged by one of these basic types plus its space dimension.
@enumx BasicType begin
    simplex = 1     # Simplicial element in any nonnegative dimension
    cube = 2        # Cube element in any nonnegative dimension
    prism = 3       # Four sided pyramid in three dimensions
    pyramid = 4     # Prism element in three dimensions
    extended = 5    # Other, more general topology, representable as topologyId
    none = 6        # Even more general topology, cannot be specified by a topologyId. Two GeometryTypes with 'none'
                    # type are equal if and only if they have the same dimension.
end

"Convert a string into a BasicType enum value."
basicType(s::String) = Dict(
    "simplex" => BasicType.simplex,
    "cube" => BasicType.cube,
    "prism" => BasicType.prism,
    "pyramid" => BasicType.pyramid,
    "extended" => BasicType.extended,
    "none" => BasicType.none)[s]

"Convert a BasicType into a string."
function Base.string(basicType::BasicType.T)
    if basicType == BasicType.simplex
        "simplex"
    elseif basicType == BasicType.cube
        "cube"
    elseif basicType == BasicType.prism
        "prism"
    elseif basicType == BasicType.pyramid
        "pyramid"
    elseif basicType == BasicType.extended
        "extended"
    else
        "none"
    end
end

"Unique label for each type of entities that can occur in a grid."
struct GeometryType
    dim::UInt8          # Dimension of the element.
    none::Bool          # Is the geometry of BasicType.none or BasicType.extended.
    topologyId::UInt32  # An identifier of the topology. If the geometry is of BasicType.none, this id must be 0.

    "Default constructor, not initializing anything."
    GeometryType() = new(0,true,0)

    "Reconstruct a Geometry type from an Id."
    GeometryType(id::UInt64) = new(UInt8(id & 0xFF), (id & 0x100) != 0, UInt32(id >> 32))

    "Constructor, using the topologyId (integer) and the dimension."
    GeometryType(topologyId::UInt32, dim::Integer) = new(dim,false,topologyId)

    "Constructor, using the topologyId (integer), the dimension and a flag for type none."
    GeometryType(topologyId::UInt32, dim::Integer, none::Bool) = new(dim,none,topologyId)

    "Map a BasicType and dimension to GeometryType"
    function GeometryType(basicType::BasicType.T, dim::Integer)
        if basicType == BasicType.simplex
            new(dim,false,0)
        elseif basicType == BasicType.cube
            new(dim,false,((dim>1) ? ((1 << dim) - 1) : 0))
        elseif basicType == BasicType.pyramid
            new(3,false,0b0011)
        elseif basicType == BasicType.prism
            new(3,false,0b0101)
        elseif basicType == BasicType.none
            new(dim,true,0)
        else
            throw(ArgumentError("GeometryType cannot be constructed from a BasicType.extended."))
        end
    end
end # GeometryType

"Construct an Id representing this GeometryType"
toId(g::GeometryType)::UInt64 = UInt64(g.dim) | (UInt64(g.none) << 8) | (UInt64(g.topologyId) << 32)

"Return the BasicType associated to the Geometry"
function basicType(g::GeometryType)
    if !g.none && (g.topologyId | 1) == 1
        BasicType.simplex
    elseif !g.none && (xor(g.topologyId, ((1 << g.dim)-1)) >> 1 == 0)
        BasicType.cube
    elseif !g.none && g.dim == 3 && (g.topologyId | 1) == 0b0101
        BasicType.prism
    elseif !g.none && g.dim == 3 && (g.topologyId | 1) == 0b0011
        BasicType.pyramid
    elseif g.none && g.topologyId == 0
        BasicType.none
    else
        BasicType.extended
    end
end

isVertex(g::GeometryType) = g.dim == 0
isLine(g::GeometryType) = g.dim == 1
isTriangle(g::GeometryType) = !g.none && g.dim == 2 && (g.topologyId | 1) == 0b0001
isQuadrilateral(g::GeometryType) = !g.none && g.dim == 2 && (g.topologyId | 1) == 0b0011
isTetrahedron(g::GeometryType) = !g.none && g.dim == 3 && (g.topologyId | 1) == 0b0001
isHexahedron(g::GeometryType) = !g.none && g.dim == 3 && (g.topologyId | 1) == 0b0111
isPyramid(g::GeometryType) = !g.none && g.dim == 3 && (g.topologyId | 1) == 0b0011
isPrism(g::GeometryType) = !g.none && g.dim == 3 && (g.topologyId | 1) == 0b0101
isSimplex(g::GeometryType) = !g.none && (g.topologyId | 1) == 1
isCube(g::GeometryType) = !g.none && (xor(g.topologyId, ((1 << g.dim)-1)) >> 1 == 0)
isNone(g::GeometryType) = g.none

"Return true if entity was constructed with a conical product in the last step."
isConical(g::GeometryType) = !g.none && (((g.topologyId & ~1) & (1 << (g.dim-1))) == 0)
"Return true if entity was constructed with a conical product in the chosen step."
isConical(g::GeometryType, step::Int) = !g.none && (((g.topologyId & ~1) & (1 << step)) == 0)
"Return true if entity was constructed with a prismatic product in the last step."
isPrismatic(g::GeometryType) = !g.none && (( (g.topologyId | 1) & (1 << (g.dim-1))) != 0)
"Return true if entity was constructed with a prismatic product in the chosen step."
isPrismatic(g::GeometryType, step::Int) = !g.none && (( (g.topologyId | 1) & (1 << step)) != 0)

"Check for equality. This method knows that in dimension 0 and 1 all BasicTypes are equal."
Base.:(==)(g1::GeometryType, g2::GeometryType) = (g1.none == g2.none) && (g1.dim == g2.dim) &&
    (g1.none || ((g1.topologyId >> 1) == (g2.topologyId >> 1)))

"Convert a GeometryType into a string."
function Base.string(g::GeometryType)
    if isVertex(g)
        "vertex"
    elseif isLine(g)
        "line"
    elseif isTriangle(g)
        "triangle"
    elseif isQuadrilateral(g)
        "quadrilateral"
    elseif isTetrahedron(g)
        "tetrahedron"
    elseif isPyramid(g)
        "pyramid"
    elseif isPrism(g)
        "prism"
    elseif isHexahedron(g)
        "hexahedron"
    else
        "none"
    end
end

"Show the GeometryType on the Standard output."
function Base.show(io::IO, g::GeometryType)
    if isSimplex(g)
        show(io, "(simplex, " * string(g.dim) * ")")
    elseif isCube(g)
        show(io, "(cube, " * string(g.dim) * ")")
    elseif isPyramid(g)
        show(io, "(pyramid, 3)")
    elseif isPrism(g)
        show(io, "(prism, 3)")
    elseif g.none && g.topologyId == 0
        show(io, "(none, " * string(g.dim) * ")")
    else
        show(io, "(other [" * bitstring(g.topologyId) * "], " * string(g.dim) * ")")
    end
end

module GeometryTypes
using DuneGeometry: BasicType,GeometryType

simplex(dim::Int) = GeometryType(BasicType.simplex, dim)
cube(dim::Int) = GeometryType(BasicType.cube, dim)
none(dim::Int) = GeometryType(BasicType.none, dim)

base(gt::GeometryType) = GeometryType(UInt32(gt.topologyId & ((1 << (gt.dim-1))-1)), gt.dim-1, gt.none)
conicalExtension(gt::GeometryType) = GeometryType(gt.topologyId, gt.dim+1, gt.none)
prismaticExtension(gt::GeometryType) = GeometryType(UInt32(gt.topologyId | ((1 << gt.dim))), gt.dim+1, gt.none)

const vertex::GeometryType = GeometryType(UInt32(0),0,false)
const line::GeometryType = GeometryType(UInt32(0),1,false)
const triangle::GeometryType = simplex(2)
const quadrilateral::GeometryType = cube(2)
const tetrahedron::GeometryType = simplex(3)
const hexahedron::GeometryType = cube(3)
const prism::GeometryType = GeometryType(BasicType.prism,3)
const pyramid::GeometryType = GeometryType(BasicType.pyramid,3)

end # module GeometryTypes