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
@test string(BasicType.simplex) == "simplex"
@test string(BasicType.cube) == "cube"
@test string(BasicType.prism) == "prism"
@test string(BasicType.pyramid) == "pyramid"
@test string(BasicType.extended) == "extended"
@test string(BasicType.none) == "none"


println("Testing GeometryType...")

vertex = GeometryType(BasicType.simplex,0)
line = GeometryType(BasicType.simplex,1)
tri = GeometryType(BasicType.simplex,2)
quad = GeometryType(BasicType.cube,2)
tet = GeometryType(BasicType.simplex,3)
hex = GeometryType(BasicType.cube,3)
prism = GeometryType(BasicType.prism,3)
pyramid = GeometryType(BasicType.pyramid,3)

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

# check isXY properties
@test isVertex(vertex)
@test isLine(line)
@test isTriangle(tri)
@test isQuadrilateral(quad)
@test isTetrahedron(tet)
@test isHexahedron(hex)
@test isPrism(prism)
@test isPyramid(pyramid)

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

# check string conversion
@test string(vertex) == "vertex"
@test string(line) == "line"
@test string(tri) == "triangle"
@test string(quad) == "quadrilateral"
@test string(tet) == "tetrahedron"
@test string(hex) == "hexahedron"
@test string(prism) == "prism"
@test string(pyramid) == "pyramid"

# check type comparison
@test vertex != line
@test vertex == GeometryType(BasicType.cube,0)
@test line == GeometryType(BasicType.cube,1)

none = GeometryType()
for gt in [vertex,line,tri,quad,tet,hex,prism,pyramid]
  @test gt == GeometryType(toId(gt)) # bug for gt == none
  if !gt.none
    @test gt == GeometryType(gt.topologyId, gt.dim)
  end
  @test gt == GeometryType(gt.topologyId, gt.dim, gt.none)
end