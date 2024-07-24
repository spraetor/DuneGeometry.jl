module Utilities

using ArgCheck
import LinearAlgebra

function AAT_L(A::AbstractArray{T,2}, ret::AbstractArray{T,2}) where {T<:Real}
  m = Base.size(A,1)
  n = Base.size(A,2)
  @argcheck Base.size(ret) == (m,m)

  for i = 1:m
    for j = 1:i
      ret[i,j] = A[i,1] * A[j,1]
      for k = 2:n
        ret[i,j] += A[i,k] * A[j,k]
      end
    end
  end
end

function cholesky_L(A::AbstractArray{T,2}, ret::AbstractArray{T,2}; checkSingular::Bool = false)::Bool where {T<:Real}
  n = Base.size(A,1)
  @argcheck Base.size(A,2) == n
  @argcheck Base.size(ret) == (n,n)

  for i = 1:n
    xDiag = copy(A[i,i])
    for j = 1:i-1
      xDiag -= ret[i,j] * ret[i,j]
    end

    # in some cases A can be singular, e.g. when checking local for
    # outside points during checkInside
    if checkSingular && ! ( xDiag > zero(T))
      return false
    end

    # otherwise this should be true always
    @assert xDiag > zero(T)
    ret[i,i] = sqrt(xDiag)

    invrii = one(T) / ret[i,i]
    for k = (i+1):n
      x = copy(A[k,i])
      for j = 1:i-1
        x -= ret[i,j] * ret[k,j]
      end
      ret[k,i] = invrii * x
    end
  end

  # return true for meaning A is non-singular
  return true
end


function invL(L::AbstractArray{T,2})::T where {T<:Real}
  n = Base.size(L,1)
  @argcheck Base.size(L,2) == n

  det = one(T)
  for i = 1:n
    det *= L[i,i];
    L[i,i] = one(T) / L[i,i]
    for j = 1:i-1
      x = L[i,j] * L[j,j]
      for k = (j+1):i-1
        x += L[i,k] * L[k,j]
      end
      L[i,j] = (-L[i,i]) * x
    end
  end
  return det
end

function LTL(L::AbstractArray{T,2}, ret::AbstractArray{T,2}) where {T<:Real}
  n = Base.size(L,1)
  @argcheck Base.size(L,2) == n
  @argcheck Base.size(ret) == (n,n)

  for i = 1:n
    for j = 1:i-1
      ret[i,j] = zero(T)
      for k = i:n
        ret[i,j] += L[k,i] * L[k,j]
      end
      ret[j,i] = ret[i,j]
    end
    ret[i,i] = zero(T)
    for k = i:n
      ret[i,i] += L[k,i] * L[k,i]
    end
  end
end

function spdInvA(A::AbstractArray{T,2})::T where {T<:Real}
  n = Base.size(A,1)
  @argcheck Base.size(A,2) == n
  L = zeros(T,n,n)
  cholesky_L(A, L)
  det = invL(L)
  LTL(L, A)
  return det
end


# Compute right pseudo-inverse of matrix `A` ant store the result in `ret`.
#
# Resturns
# The gramian determinant: `sqrt(det(A'*A))`
#
# Arguments
# - `A::AbstractArray{T,2}`: The input matrix with size `(m, n)` where `m <= n`.
# - `ref::AbstractArray{T,2}`: An output matrix with size `(n, m)`.
function rightInvA(A::AbstractArray{T,2}, ret::AbstractArray{T,2})::T where {T<:Real}
  m = Base.size(A,1)
  n = Base.size(A,2)
  @argcheck (n >= m) "Matrix has no right inverse."
  @argcheck Base.size(ret) == (n,m)
  if (n == 1) && (m == 1)
    ret[1,1] = T(1) / A[1,1]
    return abs(A[1,1])
  elseif (n == 2) && (m == 2)
    det = A[1,1] * A[2,2] - A[2,1] * A[1,2]
    detInv = T(1) / det

    ret[1,1] = A[2,2] * detInv;
    ret[2,2] = A[1,1] * detInv;
    ret[2,1] = -A[2,1] * detInv;
    ret[1,2] = -A[1,2] * detInv;
    return abs(det)
  else
    try
      ret .= pinv(A)
      if n == m
        return det(A)
      else
        return sqrt(det(A*A'))
      end
    catch e
      aat = zeros(T,m,m)
      AAT_L(A, aat)
      det = spdInvA(aat)
      ret .= A' * aat'
      return det
    end
  end
end

# Return the right pseudo-inverse of matrix A.
function rightInvA(A::AbstractArray{T,2}) where {T<:Real}
  m = Base.size(A,1)
  n = Base.size(A,2)
  @argcheck (n >= m) "Matrix has no right inverse."
  if (n == 1) && (m == 1)
    return T[one(T)/A[1,1]]
  elseif (n == 2) && (m == 2)
    det = A[1,1] * A[2,2] - A[2,1] * A[1,2]
    return T[A[2,2]/det -A[1,2]/det
            -A[2,1]/det  A[1,1]/det]
  else
    try
      return LinearAlgebra.pinv(A)
    catch e
      aat = zeros(T,m,m)
      AAT_L(A, aat)
      spdInvA(aat)
      return A' * aat'
    end
  end
end

end # module Utilities