using DuneGeometry: GeometryType, BasicType, ReferenceElement, AffineGeometry
using DuneGeometry.ReferenceElements: size, position
using Test

function isIdentity(I::AbstractArray{T,2}) where {T<:Real}
  for i=1:Base.size(I,1)
    for j=1:Base.size(I,2)
      if I[i,j] != (i==j ? T(1) : T(0))
        return false
      end
    end
  end
  return true
end

function testGeometry(geo::Geometry{T},gt::GeometryType) where {T<:Real}
  @test affine(geo)
  @test AffineGeometries.type(geo) == gt

  ref = referenceElement(geo)
  @test ReferenceElements.type(ref) == gt
  @test geo.mydimension == gt.dim

  @test corners(geo) == size(ref,geo.mydimension)
  @test center(geo) == position(ref,1,0)

  for i=1:corners(geo)
    @test corner(geo,i) == position(ref,i,geo.mydimension)
  end

  # test specific geometries
  if isLine(gt)
    @test corner(geo, 1) == T[0.0]
    @test corner(geo, 2) == T[1.0]

  elseif isTriangle(gt)
    @test corner(geo, 1) == T[0.0,0.0]
    @test corner(geo, 2) == T[1.0,0.0]
    @test corner(geo, 3) == T[0.0,1.0]

  elseif isQuadrilateral(gt)
    @test corner(geo, 1) == T[0.0,0.0]
    @test corner(geo, 2) == T[1.0,0.0]
    @test corner(geo, 3) == T[0.0,1.0]
    @test corner(geo, 4) == T[1.0,1.0]
  end

  points = geo.mydimension == 1 ? [ T[i/6.0] for i=0:5 ] :
           geo.mydimension == 2 ? [ T[i/6.0,j/6.0] for i=0:5 for j=0:5 ] :
           geo.mydimension == 3 ? [ T[i/6.0,j/6.0,k/6.0] for i=0:5 for j=0:5 for k=0:5 ] :
                                  Vector{Vector{T}}(undef,0)

  points = points[ ReferenceElements.checkInside.(Ref(ref), points) ]
  for p in points
    @test toGlobal(geo,p) == p
    @test toLocal(geo,p) == p
  end

  c = position(ref,1,0)
  @test AffineGeometries.volume(geo) == integrationElement(geo, c) * ReferenceElements.volume(referenceElement(geo))
  @test isIdentity(jacobianInverseTransposed(geo,c) * jacobianTransposed(geo,c))
  @test isIdentity(jacobian(geo,c) * jacobianInverse(geo,c))
end



geometryTypes = [GeometryType(BasicType.simplex,1),
                 GeometryType(BasicType.simplex,2), GeometryType(BasicType.cube,2),
                 GeometryType(BasicType.simplex,3), GeometryType(BasicType.cube,3),
                 GeometryType(BasicType.prism,3), GeometryType(BasicType.pyramid,3)]

for gt in geometryTypes
  for T in (Float32,Float64,BigFloat)
    println("Testing AffineGeometry{$(T)} $(gt)")
    ref = ReferenceElement{T}(gt)
    testGeometry(ReferenceElements.geometry(AffineGeometry{T},ref,1,0),gt)

    # test constructions
    if isTriangle(gt)
      origin = position(ref,1,2)
      jt = [position(ref,2,2)-origin position(ref,3,2)-origin]
      geo1 = AffineGeometry{T}(ref, [position(ref,i,2) for i=1:size(ref,2)])
      geo2 = AffineGeometry{T}(tri, [position(ref,i,2) for i=1:size(ref,2)])
      geo3 = AffineGeometry{T}(ref, origin, jt)
      geo4 = AffineGeometry{T}(tri, origin, jt)

      testGeometry(geo1,gt)
      testGeometry(geo3,gt)
    end

    # @test geo1 == geo2
    # @test get1 == geo3
    # @test get1 == geo4
  end
end