using DuneGeometry: GeometryType, BasicType, ReferenceElement, AffineGeometry
using DuneGeometry: size, position

geometryTypes = [GeometryType(BasicType.simplex,1),
                 GeometryType(BasicType.simplex,2), GeometryType(BasicType.cube,2),
                 GeometryType(BasicType.simplex,3), GeometryType(BasicType.cube,3),
                 GeometryType(BasicType.prism,3), GeometryType(BasicType.pyramid,3)]

for gt in geometryTypes
  for T in (Float32,Float64,BigFloat)
    println("Testing AffineGeometry{$(T)} $(gt)")
    ref = ReferenceElement{T}(gt)
    testGeometry(geometry(AffineGeometry{T},ref,1,0),gt)

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
  end
end