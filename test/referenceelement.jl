using DuneGeometry: GeometryType, BasicType, ReferenceElement
using Test

import DuneGeometry.ReferenceElements

println("Testing ReferenceElement")

tri = GeometryType(BasicType.simplex,2)
ref = ReferenceElement{Float64}(tri)

@test ReferenceElements.size(ref,0) == 1
@test ReferenceElements.size(ref,1) == 3
@test ReferenceElements.size(ref,2) == 3

# subentities of face
@test ReferenceElements.size(ref,1,0,0) == 1
@test ReferenceElements.size(ref,1,0,1) == 3
@test ReferenceElements.size(ref,1,0,2) == 3

# subentities of edge
@test ReferenceElements.size(ref,1,1,1) == 1
@test ReferenceElements.size(ref,1,1,2) == 2
@test ReferenceElements.size(ref,2,1,1) == 1
@test ReferenceElements.size(ref,2,1,2) == 2
@test ReferenceElements.size(ref,3,1,1) == 1
@test ReferenceElements.size(ref,3,1,2) == 2

# subentities of vertex
@test ReferenceElements.size(ref,1,2,2) == 1
@test ReferenceElements.size(ref,2,2,2) == 1
@test ReferenceElements.size(ref,3,2,2) == 1
