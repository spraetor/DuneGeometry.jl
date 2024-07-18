using DuneGeometry: GeometryType, BasicType
using Test

println("Testing ReferenceElementsImpl")
import DuneGeometry.ReferenceElementsImpl

# 1. check size and subTopologyId

vertex = GeometryType(BasicType.simplex,0)
@test ReferenceElementsImpl.size(vertex.topologyId, vertex.dim, 0) == 1
@test GeometryType(ReferenceElementsImpl.subTopologyId(vertex.topologyId, vertex.dim, 0, 1), 0) == vertex

line = GeometryType(BasicType.simplex,1)
@test ReferenceElementsImpl.size(line.topologyId, line.dim, 0) == 1
@test ReferenceElementsImpl.size(line.topologyId, line.dim, 1) == 2
@test GeometryType(ReferenceElementsImpl.subTopologyId(line.topologyId, line.dim, 0, 1), 1) == line
@test GeometryType(ReferenceElementsImpl.subTopologyId(line.topologyId, line.dim, 1, 1), 0) == vertex
@test GeometryType(ReferenceElementsImpl.subTopologyId(line.topologyId, line.dim, 1, 2), 0) == vertex

tri = GeometryType(BasicType.simplex,2)
@test ReferenceElementsImpl.size(tri.topologyId, tri.dim, 0) == 1
@test ReferenceElementsImpl.size(tri.topologyId, tri.dim, 1) == 3
@test ReferenceElementsImpl.size(tri.topologyId, tri.dim, 2) == 3
@test GeometryType(ReferenceElementsImpl.subTopologyId(tri.topologyId, tri.dim, 0, 1), 2) == tri
@test GeometryType(ReferenceElementsImpl.subTopologyId(tri.topologyId, tri.dim, 1, 1), 1) == line
@test GeometryType(ReferenceElementsImpl.subTopologyId(tri.topologyId, tri.dim, 1, 2), 1) == line
@test GeometryType(ReferenceElementsImpl.subTopologyId(tri.topologyId, tri.dim, 1, 3), 1) == line
@test GeometryType(ReferenceElementsImpl.subTopologyId(tri.topologyId, tri.dim, 2, 1), 0) == vertex
@test GeometryType(ReferenceElementsImpl.subTopologyId(tri.topologyId, tri.dim, 2, 2), 0) == vertex
@test GeometryType(ReferenceElementsImpl.subTopologyId(tri.topologyId, tri.dim, 2, 3), 0) == vertex

tet = GeometryType(BasicType.simplex,3)
@test ReferenceElementsImpl.size(tet.topologyId, tet.dim, 0) == 1
@test ReferenceElementsImpl.size(tet.topologyId, tet.dim, 1) == 4
@test ReferenceElementsImpl.size(tet.topologyId, tet.dim, 2) == 6
@test ReferenceElementsImpl.size(tet.topologyId, tet.dim, 3) == 4

quad = GeometryType(BasicType.cube,2)
@test ReferenceElementsImpl.size(quad.topologyId, quad.dim, 0) == 1
@test ReferenceElementsImpl.size(quad.topologyId, quad.dim, 1) == 4
@test ReferenceElementsImpl.size(quad.topologyId, quad.dim, 2) == 4

hex = GeometryType(BasicType.cube,3)
@test ReferenceElementsImpl.size(hex.topologyId, hex.dim, 0) == 1
@test ReferenceElementsImpl.size(hex.topologyId, hex.dim, 1) == 6
@test ReferenceElementsImpl.size(hex.topologyId, hex.dim, 2) == 12
@test ReferenceElementsImpl.size(hex.topologyId, hex.dim, 3) == 8

prism = GeometryType(BasicType.prism,3)
@test ReferenceElementsImpl.size(prism.topologyId, prism.dim, 0) == 1
@test ReferenceElementsImpl.size(prism.topologyId, prism.dim, 1) == 5
@test ReferenceElementsImpl.size(prism.topologyId, prism.dim, 2) == 9
@test ReferenceElementsImpl.size(prism.topologyId, prism.dim, 3) == 6

pyramid = GeometryType(BasicType.pyramid,3)
@test ReferenceElementsImpl.size(pyramid.topologyId, pyramid.dim, 0) == 1
@test ReferenceElementsImpl.size(pyramid.topologyId, pyramid.dim, 1) == 5
@test ReferenceElementsImpl.size(pyramid.topologyId, pyramid.dim, 2) == 8
@test ReferenceElementsImpl.size(pyramid.topologyId, pyramid.dim, 3) == 5

# 2. checkInside

@test ReferenceElementsImpl.checkInside(tri.topologyId, tri.dim, Float64[0.1,0.1], eps(Float64))
@test ReferenceElementsImpl.checkInside(tri.topologyId, tri.dim, Float64[0.0,0.0], eps(Float64))
@test ReferenceElementsImpl.checkInside(tri.topologyId, tri.dim, Float64[1.0,0.0], eps(Float64))
@test ReferenceElementsImpl.checkInside(tri.topologyId, tri.dim, Float64[0.0,1.0], eps(Float64))
@test !ReferenceElementsImpl.checkInside(tri.topologyId, tri.dim, Float64[-0.1,0.1], eps(Float64))

# 3. check referenceCorners!

corners_vertex = [Float64[]]
ReferenceElementsImpl.referenceCorners!(vertex.topologyId, vertex.dim, corners_vertex)
@test Float64[] in corners_vertex

corners_line = [ Vector{Float64}(undef,1) for _ in 1:2 ]
ReferenceElementsImpl.referenceCorners!(line.topologyId, line.dim, corners_line)
@test Float64[0] in corners_line
@test Float64[1] in corners_line

corners_tri = [ Vector{Float64}(undef,2) for _ in 1:3 ]
ReferenceElementsImpl.referenceCorners!(tri.topologyId, tri.dim, corners_tri)
@test Float64[0,0] in corners_tri
@test Float64[1,0] in corners_tri
@test Float64[0,1] in corners_tri

corners_tet = [ Vector{Float64}(undef,3) for _ in 1:4 ]
ReferenceElementsImpl.referenceCorners!(tet.topologyId, tet.dim, corners_tet)
@test Float64[0,0,0] in corners_tet
@test Float64[1,0,0] in corners_tet
@test Float64[0,1,0] in corners_tet
@test Float64[0,0,1] in corners_tet

corners_quad = [ Vector{Float64}(undef,2) for _ in 1:4 ]
ReferenceElementsImpl.referenceCorners!(quad.topologyId, quad.dim, corners_quad)
@test Float64[0,0] in corners_quad
@test Float64[1,0] in corners_quad
@test Float64[0,1] in corners_quad
@test Float64[1,1] in corners_quad

corners_hex = [ Vector{Float64}(undef,3) for _ in 1:8 ]
ReferenceElementsImpl.referenceCorners!(hex.topologyId, hex.dim, corners_hex)
@test Float64[0,0,0] in corners_hex
@test Float64[1,0,0] in corners_hex
@test Float64[0,1,0] in corners_hex
@test Float64[1,1,0] in corners_hex
@test Float64[0,0,1] in corners_hex
@test Float64[1,0,1] in corners_hex
@test Float64[0,1,1] in corners_hex
@test Float64[1,1,1] in corners_hex

corners_prism = [ Vector{Float64}(undef,3) for _ in 1:6 ]
ReferenceElementsImpl.referenceCorners!(prism.topologyId, prism.dim, corners_prism)
@test Float64[0,0,0] in corners_prism
@test Float64[1,0,0] in corners_prism
@test Float64[0,1,0] in corners_prism
@test Float64[0,0,1] in corners_prism
@test Float64[1,0,1] in corners_prism
@test Float64[0,1,1] in corners_prism

corners_pyramid = [ Vector{Float64}(undef,3) for _ in 1:5 ]
ReferenceElementsImpl.referenceCorners!(pyramid.topologyId, pyramid.dim, corners_pyramid)
@test Float64[0,0,0] in corners_pyramid
@test Float64[1,0,0] in corners_pyramid
@test Float64[0,1,0] in corners_pyramid
@test Float64[1,1,0] in corners_pyramid
@test Float64[0,0,1] in corners_pyramid

# 3. check volume

@test ReferenceElementsImpl.referenceVolume(Float64, vertex.topologyId, vertex.dim) == 1.0
@test ReferenceElementsImpl.referenceVolume(Float64, line.topologyId, line.dim) == 1.0
@test ReferenceElementsImpl.referenceVolume(Float64, tri.topologyId, tri.dim) == 0.5
@test ReferenceElementsImpl.referenceVolume(Float64, tet.topologyId, tet.dim) == 1.0/6.0
@test ReferenceElementsImpl.referenceVolume(Float64, quad.topologyId, quad.dim) == 1.0
@test ReferenceElementsImpl.referenceVolume(Float64, hex.topologyId, hex.dim) == 1.0
@test ReferenceElementsImpl.referenceVolume(Float64, prism.topologyId, prism.dim) == 0.5
@test ReferenceElementsImpl.referenceVolume(Float64, pyramid.topologyId, pyramid.dim) == 1.0/3.0

# 4. check referenceOrigins!

origins_tri0 = [ Vector{Float64}(undef,2) for _ in 1:1 ]
origins_tri1 = [ Vector{Float64}(undef,2) for _ in 1:3 ]
origins_tri2 = [ Vector{Float64}(undef,2) for _ in 1:3 ]
ReferenceElementsImpl.referenceOrigins!(tri.topologyId, tri.dim, 0, origins_tri0)
ReferenceElementsImpl.referenceOrigins!(tri.topologyId, tri.dim, 1, origins_tri1)
ReferenceElementsImpl.referenceOrigins!(tri.topologyId, tri.dim, 2, origins_tri2)
@test Float64[0,0] in origins_tri0
@test Float64[0,0] in origins_tri1
@test Float64[0,0] in origins_tri2

origins_tet0 = [ Vector{Float64}(undef,3) for _ in 1:1 ]
origins_tet1 = [ Vector{Float64}(undef,3) for _ in 1:4 ]
origins_tet2 = [ Vector{Float64}(undef,3) for _ in 1:6 ]
origins_tet3 = [ Vector{Float64}(undef,3) for _ in 1:4 ]
ReferenceElementsImpl.referenceOrigins!(tet.topologyId, tet.dim, 0, origins_tet0)
ReferenceElementsImpl.referenceOrigins!(tet.topologyId, tet.dim, 1, origins_tet1)
ReferenceElementsImpl.referenceOrigins!(tet.topologyId, tet.dim, 2, origins_tet2)
ReferenceElementsImpl.referenceOrigins!(tet.topologyId, tet.dim, 3, origins_tet3)
@test Float64[0,0,0] in origins_tet0
@test Float64[0,0,0] in origins_tet1
@test Float64[0,0,0] in origins_tet2
@test Float64[0,0,0] in origins_tet3

origins_prism0 = [ Vector{Float64}(undef,3) for _ in 1:1 ]
origins_prism1 = [ Vector{Float64}(undef,3) for _ in 1:5 ]
origins_prism2 = [ Vector{Float64}(undef,3) for _ in 1:9 ]
origins_prism3 = [ Vector{Float64}(undef,3) for _ in 1:6 ]
ReferenceElementsImpl.referenceOrigins!(prism.topologyId, prism.dim, 0, origins_prism0)
ReferenceElementsImpl.referenceOrigins!(prism.topologyId, prism.dim, 1, origins_prism1)
ReferenceElementsImpl.referenceOrigins!(prism.topologyId, prism.dim, 2, origins_prism2)
ReferenceElementsImpl.referenceOrigins!(prism.topologyId, prism.dim, 3, origins_prism3)
@test Float64[0,0,0] in origins_prism0
@test Float64[0,0,0] in origins_prism1
@test Float64[0,0,0] in origins_prism2
@test Float64[0,0,0] in origins_prism3


# 5. check referenceEmbeddings!


origins_tri0_ = [ Vector{Float64}(undef,2) for _ in 1:1 ]
origins_tri1_ = [ Vector{Float64}(undef,2) for _ in 1:3 ]
origins_tri2_ = [ Vector{Float64}(undef,2) for _ in 1:3 ]
jacobianTransposeds_tri0 = [ Matrix{Float64}(undef,2,2) for _ in 1:1 ]
jacobianTransposeds_tri1 = [ Matrix{Float64}(undef,2,2) for _ in 1:3 ]
jacobianTransposeds_tri2 = [ Matrix{Float64}(undef,2,2) for _ in 1:3 ]
ReferenceElementsImpl.referenceEmbeddings!(tri.topologyId, tri.dim, 0, origins_tri0_, jacobianTransposeds_tri0)
ReferenceElementsImpl.referenceEmbeddings!(tri.topologyId, tri.dim, 1, origins_tri1_, jacobianTransposeds_tri1)
ReferenceElementsImpl.referenceEmbeddings!(tri.topologyId, tri.dim, 2, origins_tri2_, jacobianTransposeds_tri2)
@test origins_tri0_ == origins_tri0
@test origins_tri1_ == origins_tri1
@test origins_tri2_ == origins_tri2

for jt in jacobianTransposeds_tri0
  @test jt == Float64[1 0; 0 1]
end

for jt in jacobianTransposeds_tri2
  @test jt == zeros(Float64, 2,2)
end
