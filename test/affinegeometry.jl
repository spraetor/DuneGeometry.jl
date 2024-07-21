using DuneGeometry: GeometryType, BasicType, ReferenceElement, AffineGeometry
using DuneGeometry.ReferenceElements: size, position
using Test

println("Testing AffineGeometry")

let tri = GeometryType(BasicType.simplex,2), ref = ReferenceElement{Float64}(tri)

  geo = AffineGeometry{Float64}(ref, [position(ref,i,2) for i=1:size(ref,2)])
  let
    origin = position(ref,1,2)
    jt = [position(ref,2,2)-origin position(ref,3,2)-origin]
    geo2 = AffineGeometry{Float64}(tri, [position(ref,i,2) for i=1:size(ref,2)])
    geo3 = AffineGeometry{Float64}(ref, origin, jt)
    geo4 = AffineGeometry{Float64}(tri, origin, jt)
  end

  # @test geo1 == geo2
  # @test get1 == geo3
  # @test get1 == geo4

  @test affine(geo)
  @test AffineGeometries.type(geo) == tri
  @test corners(geo) == 3

  @test corner(geo, 1) == Float64[0.0,0.0]
  @test corner(geo, 2) == Float64[1.0,0.0]
  @test corner(geo, 3) == Float64[0.0,1.0]

  @test center(geo) == Float64[1.0,1.0] ./ 3.0

  points = [ Float64[i/6.0, j/6.0] for i=0:5 for j=0:5 ]
  points = points[ [p[1] + p[2] <= 1.0 for p in points] ]
  for p in points
    @test toGlobal(geo,p) == p
    @test toLocal(geo,p) == p
  end

end