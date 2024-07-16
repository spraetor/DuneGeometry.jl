using QuadRules: GeometryType, BasicType, ReferenceElement
using Test

import QuadRules.ReferenceElements

println("Testing ReferenceElement")

tri = GeometryType(BasicType.simplex,2)
ref = ReferenceElement{Float64}(tri)

@test ReferenceElements.size(ref,0) == 1
@test ReferenceElements.size(ref,1) == 3
@test ReferenceElements.size(ref,2) == 3