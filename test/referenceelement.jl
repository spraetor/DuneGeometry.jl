using DuneGeometry: GeometryType, BasicType, ReferenceElement
using DuneGeometry.ReferenceElements: size, position, type, volume
using DuneGeometry.AffineGeometries: AffineGeometry,affine,integrationElement
using Test

function testTriangleReferenceElement(ref::ReferenceElement{T}, tri::GeometryType) where {T<:Real}
  @test size(ref,0) == 1
  @test size(ref,1) == 3
  @test size(ref,2) == 3

  # subentities of face
  @test size(ref,1,0,0) == 1
  @test size(ref,1,0,1) == 3
  @test size(ref,1,0,2) == 3

  # subentities of edge
  @test size(ref,1,1,1) == 1
  @test size(ref,1,1,2) == 2
  @test size(ref,2,1,1) == 1
  @test size(ref,2,1,2) == 2
  @test size(ref,3,1,1) == 1
  @test size(ref,3,1,2) == 2

  # subentities of vertex
  @test size(ref,1,2,2) == 1
  @test size(ref,2,2,2) == 1
  @test size(ref,3,2,2) == 1


  # test sub entities
  @test subEntity(ref, 1,0,1,0) == 1

  for cc = 0:2
    for ii = 1:size(ref,1,0,cc)
      @test subEntity(ref,1,0,ii,cc) == ii
    end
  end

  for c = 0:ref.dimension
    for i = 1:size(ref,c)
      for cc = c:ref.dimension
        subEntityRange,contains = subEntities(ref,i,c,cc)
        @test length(subEntityRange) == size(ref,i,c,cc)
        for ii in eachindex(subEntityRange)
          @test subEntityRange[ii] == subEntity(ref,i,c,ii,cc)
          @test contains[subEntityRange[ii]]
        end
      end
    end
  end


  # test geometry types
  @test type(ref,1,0) == type(ref)
  @test type(ref,1,0) == tri
  for i = 1:size(ref,1)
    @test type(ref,i,1) == GeometryType(BasicType.simplex,1)
  end
  for i = 1:size(ref,2)
    @test type(ref,i,2) == GeometryType(BasicType.simplex,0)
  end


  # test positions
  @test position(ref,1,0) == Float64[1.0/3.0, 1.0/3.0]
  @test position(ref,1,1) == Float64[0.5,0.0]
  @test position(ref,2,1) == Float64[0.0,0.5]
  @test position(ref,3,1) == Float64[0.5,0.5]
  @test position(ref,1,2) == Float64[0.0,0.0]
  @test position(ref,2,2) == Float64[1.0,0.0]
  @test position(ref,3,2) == Float64[0.0,1.0]


  # test checkInside
  @test checkInside(ref, position(ref,1,0))
  @test checkInside(ref, Float64[0.1, 0.1])
  @test checkInside(ref, Float64[0.0, 0.0])
  @test checkInside(ref, Float64[1.0, 0.0])
  @test checkInside(ref, Float64[0.0, 1.0])
  @test !checkInside(ref, Float64[-1.0, -1.0])


  # test geometry
  geo0 = geometry(AffineGeometry{Float64}, ref, 1, 0)
  geo1 = [geometry(AffineGeometry{Float64}, ref, i, 1) for i = 1:3]
  geo2 = [geometry(AffineGeometry{Float64}, ref, i, 2) for i = 1:3]

  @test affine(geo0)
  for i = 1:3
    @test affine(geo1[i])
  end
  for i = 1:3
    @test affine(geo2[i])
  end


  # test volume
  @test volume(ref) == Float64(0.5)


  # test integrationOuterNormal
  @test integrationOuterNormal(ref,1) == Float64[0.0, -1.0]
  @test integrationOuterNormal(ref,2) == Float64[-1.0, 0.0]
  @test integrationOuterNormal(ref,3) == Float64[1.0, 1.0]

  for i = 1:3
    n = integrationOuterNormal(ref,i)
    @test sqrt(sum(n.*n)) == integrationElement(geo1[i], position(ref,i,1))
  end
end



# import DuneGeometry.ReferenceElements
let tri = GeometryType(BasicType.simplex,2)
  for T in (Float32,Float64,BigFloat)
    println("Testing ReferenceElement{$(T)}")

    ref = ReferenceElement{Float64}(tri)
    testTriangleReferenceElement(ref,tri)
  end
end
