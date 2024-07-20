using DuneGeometry: GeometryType, BasicType, ReferenceElement
using DuneGeometry.ReferenceElements: size, position
using Test

# import DuneGeometry.ReferenceElements

println("Testing ReferenceElement")

tri = GeometryType(BasicType.simplex,2)
ref = ReferenceElement{Float64}(tri)

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
# @test position(ref,1,1) == Float64[0.5,0.0]
# @test position(ref,2,1) == Float64[0.0,0.5]
# @test position(ref,3,1) == Float64[0.5,0.5]
@test position(ref,1,2) == Float64[0.0,0.0]
@test position(ref,2,2) == Float64[1.0,0.0]
@test position(ref,3,2) == Float64[0.0,1.0]

# NOTE: Something with SubEntityInfo.numbering is still wrong.