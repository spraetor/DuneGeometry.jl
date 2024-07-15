using QuadRules: GeometryType, BasicType
using Test

println("Testing ReferenceElementsImpl")
import QuadRules.ReferenceElementsImpl

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
