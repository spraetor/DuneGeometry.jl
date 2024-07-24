using Test
import DuneGeometry.Utilities
import LinearAlgebra

println("Testing Utilities")

function testUtilities(A::AbstractArray{T,2}) where {T<:Real}
  AAt = A * A'

  # test AAT_L
  AAt_L = zeros(eltype(AAt), Base.size(AAt))
  Utilities.AAT_L(A, AAt_L)

  for i in axes(AAt,1)
    for j = 1:i
      @test AAt[i,j] ≈ AAt_L[i,j]
    end
  end

  # test cholesky_L
  C = LinearAlgebra.cholesky(AAt)
  C_L = zeros(eltype(AAt), Base.size(AAt))
  issingular = Utilities.cholesky_L(AAt,C_L; checkSingular=true)

  @test LinearAlgebra.issuccess(C) == issingular
  for i in axes(AAt,1)
    for j = 1:i
      @test C.L[i,j] ≈ C_L[i,j]
    end
  end

  # test invL
  AAt_inv = inv(AAt)
  det_AAt = LinearAlgebra.det(AAt)
  C_L_inv = copy(C_L)
  det_C_L = Utilities.invL(C_L_inv)
  AAt_L_inv = C_L_inv' * C_L_inv

  @test sqrt(det_AAt) ≈ det_C_L
  for i in axes(AAt_inv,1)
    for j = 1:i
      @test AAt_inv[i,j] ≈ AAt_L_inv[i,j]
    end
  end

  # test spdInvA
  AAt_L2 = copy(AAt_L)
  det_AAt_L2 = Utilities.spdInvA(AAt_L2)

  @test sqrt(det_AAt) ≈ det_AAt_L2
  for i in axes(AAt_inv,1)
    for j in axes(AAt_inv,2)
      @test AAt_inv[i,j] ≈ AAt_L2[i,j]
    end
  end

  # test rightInvA
  A_pinv = LinearAlgebra.pinv(A)
  A_rightinv = zeros(eltype(A_pinv), Base.size(A_pinv))
  sqrt_det_AAt = Utilities.rightInvA(A, A_rightinv)

  @test sqrt(det_AAt) ≈ sqrt_det_AAt
  for i in axes(A_pinv,1)
    for j in axes(A_pinv,2)
      @test A_pinv[i,j] ≈ A_rightinv[i,j]
    end
  end
end



A1 = Float64[1.0;;]
A2 = Float64[1 2; 4 0.2]
A3 = Float64[1 2 3; 4 0.2 6]

for A in (A1,A2,A3)
  testUtilities(A)
end