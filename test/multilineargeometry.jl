using DuneGeometry: GeometryType, BasicType, ReferenceElement, MultiLinearGeometry
using DuneGeometry: size, position

geometryTypes = [GeometryType(BasicType.simplex,1),
                 GeometryType(BasicType.simplex,2), GeometryType(BasicType.cube,2),
                 GeometryType(BasicType.simplex,3), GeometryType(BasicType.cube,3),
                 GeometryType(BasicType.prism,3), GeometryType(BasicType.pyramid,3)]

for gt in geometryTypes
  for T in (Float32,Float64,BigFloat)
    println("Testing MultiLinearGeometry{$(T)} $(gt)")
    ref = ReferenceElement{T}(gt)

    testGeometry(MultiLinearGeometry{T}(ref,[position(ref,i,gt.dim) for i = 1:size(ref,gt.dim)]),gt)
  end
end