using DuneGeometry: GeometryType, BasicType, ReferenceElement, AffineGeometry
using DuneGeometry.ReferenceElements: size, position
using Test

println("Testing AffineGeometry")

let tri = GeometryType(BasicType.simplex,2), ref = ReferenceElement{Float64}(tri)

  geo1 = AffineGeometry{Float64}(ref, [position(ref,i,2) for i=1:size(ref,2)])
  geo2 = AffineGeometry{Float64}(tri, [position(ref,i,2) for i=1:size(ref,2)])

  origin = position(ref,1,2)
  jt = [position(ref,2,2)-origin position(ref,3,2)-origin]
  geo3 = AffineGeometry{Float64}(ref, origin, jt)
  geo4 = AffineGeometry{Float64}(tri, origin, jt)

  # @test geo1 == geo2
  # @test get1 == geo3
  # @test get1 == geo4

  @test affine(geo1)
  @test AffineGeometries.type(geo1) == tri
  @test corners(geo1) == 3
end