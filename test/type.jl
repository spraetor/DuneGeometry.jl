using DuneGeometry
using Test

println("Testing BasicType")

# check conversion of string to BasicType
@test basicType("simplex") == BasicType.simplex
@test basicType("cube") == BasicType.cube
@test basicType("prism") == BasicType.prism
@test basicType("pyramid") == BasicType.pyramid
@test basicType("extended") == BasicType.extended
@test basicType("none") == BasicType.none

# check convertion of BAsicType to string
@test toString(BasicType.simplex) == "simplex"
@test toString(BasicType.cube) == "cube"
@test toString(BasicType.prism) == "prism"
@test toString(BasicType.pyramid) == "pyramid"
@test toString(BasicType.extended) == "extended"
@test toString(BasicType.none) == "none"


println("Testing GeometryType...")

vertex = GeometryType(BasicType.simplex,0)
line = GeometryType(BasicType.simplex,1)
tri = GeometryType(BasicType.simplex,2)
quad = GeometryType(BasicType.cube,2)
tet = GeometryType(BasicType.simplex,3)
hex = GeometryType(BasicType.cube,3)
prism = GeometryType(BasicType.prism,3)
pyramid = GeometryType(BasicType.pyramid,3)
none0 = GeometryType(BasicType.none,0)
none1 = GeometryType(BasicType.none,1)
none2 = GeometryType(BasicType.none,2)

let correctException = false
  try
    ext = GeometryType(BasicType.extended,1)
  catch ArgumentError
    correctException = true
  end
  @test correctException
end

show(devnull, none0)
show(devnull, none1)
show(devnull, none2)

default = GeometryType()
@test default == none0
@test default != none1

# check conversion to BasicType
@test basicType(vertex) == BasicType.simplex
@test basicType(line) == BasicType.simplex
@test basicType(tri) == BasicType.simplex
@test basicType(tet) == BasicType.simplex
@test basicType(quad) == BasicType.cube
@test basicType(hex) == BasicType.cube
@test basicType(quad) == BasicType.cube
@test basicType(prism) == BasicType.prism
@test basicType(pyramid) == BasicType.pyramid
@test basicType(none0) == BasicType.none
@test basicType(none1) == BasicType.none
@test basicType(none2) == BasicType.none

# check isXY properties
@test isVertex(vertex)
@test isLine(line)
@test isTriangle(tri)
@test isQuadrilateral(quad)
@test isTetrahedron(tet)
@test isHexahedron(hex)
@test isPrism(prism)
@test isPyramid(pyramid)
@test isNone(none0)
@test isNone(none1)
@test isNone(none2)
@test !isNone(vertex)

@test isSimplex(vertex)
@test isSimplex(line)
@test isSimplex(tri)
@test isSimplex(tet)

@test isCube(vertex)
@test isCube(line)
@test isCube(quad)
@test isCube(hex)

@test isPrismatic(line)
@test isPrismatic(quad)
@test isPrismatic(hex)
@test isPrismatic(prism)

@test isConical(line)
@test isConical(tri)
@test isConical(tet)
@test isConical(pyramid)

@test isPrismatic(quad,1)
@test isPrismatic(hex,2)
@test isConical(tri,1)
@test isConical(tet,2)

# check string conversion
@test toString(vertex) == "vertex"
@test toString(line) == "line"
@test toString(tri) == "triangle"
@test toString(quad) == "quadrilateral"
@test toString(tet) == "tetrahedron"
@test toString(hex) == "hexahedron"
@test toString(prism) == "prism"
@test toString(pyramid) == "pyramid"

# check type comparison
@test vertex != line
@test vertex == GeometryType(BasicType.cube,0)
@test line == GeometryType(BasicType.cube,1)

@test none0 != none1
@test none0 != none2
@test none1 != none2

# check GeometryType names
import DuneGeometry.GeometryTypes
@test vertex == GeometryTypes.vertex
@test line == GeometryTypes.line
@test tri == GeometryTypes.triangle
@test tet == GeometryTypes.tetrahedron
@test quad == GeometryTypes.quadrilateral
@test hex == GeometryTypes.hexahedron
@test prism == GeometryTypes.prism
@test pyramid == GeometryTypes.pyramid
@test none0 == GeometryTypes.none(0)
@test none1 == GeometryTypes.none(1)
@test none2 == GeometryTypes.none(2)

import DuneGeometry.TypesImpl
@test GeometryTypes.base(tet) == GeometryType(TypesImpl.baseTopologyId(tet.topologyId,tet.dim,0), tet.dim-1, false)
@test GeometryTypes.conicalExtension(line) == tri
@test GeometryTypes.conicalExtension(tri) == tet
@test GeometryTypes.conicalExtension(quad) == pyramid
@test GeometryTypes.prismaticExtension(line) == quad
@test GeometryTypes.prismaticExtension(tri) == prism
@test GeometryTypes.prismaticExtension(quad) == hex

for gt in [vertex,line,tri,quad,tet,hex,prism,pyramid,none0,none1,none2]
  @test gt == GeometryType(toId(gt))
  if !gt.none
    @test gt == GeometryType(gt.topologyId, gt.dim)
  end
  @test gt == GeometryType(gt.topologyId, gt.dim, gt.none)
end