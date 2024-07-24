using DuneGeometry: AbstractGeometry, GeometryType, size, position
using Test

function isIdentity(I::AbstractArray{T,2}) where {T<:Real}
  for i=1:Base.size(I,1)
    for j=1:Base.size(I,2)
      if (I[i,j] - (i==j ? T(1) : T(0))) > 16*eps(T)
        return false
      end
    end
  end
  return true
end

function testGeometry(geo::AbstractGeometry{T},gt::GeometryType) where {T<:Real}
  @test type(geo) == gt

  ref = referenceElement(geo)
  @test type(ref) == gt
  @test geo.mydimension == gt.dim

  @test corners(geo) == size(ref,geo.mydimension)
  @test maximum(abs.(center(geo) .- position(ref,1,0))) < 16*eps(T)

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

  points = points[ checkInside.(Ref(ref), points) ]
  for p in points
    @test maximum(abs.(toLocal(geo,toGlobal(geo,p)) .- p)) < 16*eps(T)
  end

  c = position(ref,1,0)
  @test volume(geo) == integrationElement(geo, c) * volume(referenceElement(geo))

  @test isIdentity(jacobianInverseTransposed(geo,c) * jacobianTransposed(geo,c))
  @test isIdentity(jacobian(geo,c) * jacobianInverse(geo,c))
end
