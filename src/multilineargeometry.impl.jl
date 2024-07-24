module MultiLinearGeometryImpl

import ..TypesImpl
using ArgCheck

mutable struct CornerStorage{T<:Real}
  corners::Vector{Vector{T}}
  state::Int

  CornerStorage(c::Vector{Vector{T}}, s::Integer) where {T<:Real} = new{T}(c, Int(s))
  CornerStorage(c::Vector{Vector{T}}) where {T<:Real} = new{T}(c, 1)
end

function get(cs::CornerStorage{T}) where {T<:Real}
  cs.corners[cs.state]
end


function toGlobal(add::Bool, topologyId::UInt32, mydim::Integer, dim::Integer, cs::CornerStorage{T}, df::T,
                  x::AbstractVector{T}, rf::T, X::AbstractVector{T}) where {T<:Real}
  @argcheck length(cs.corners) > 0
  @argcheck length(X) == length(cs.corners[1])
  if dim == 0
    corner = get(cs)
    for i in eachindex(X)
      X[i] = (add ? X[i] + rf*corner[i] : rf*corner[i])
    end
    cs.state += 1
  else
    tol = 16 * eps(T)
    xn = df * x[dim]
    cxn = T(1) - xn

    if TypesImpl.isPrism(topologyId, mydim, mydim-dim)
      # apply (1-xn) times mapping for bottom
      toGlobal(add, topologyId, mydim, dim-1, cs, df, x, rf*cxn, X)
      # apply xn times mapping for top
      toGlobal(true, topologyId, mydim, dim-1, cs, df, x, rf*xn, X)
    else
      @assert TypesImpl.isPyramid(topologyId, mydim, mydim-dim)
      # apply (1-xn) times mapping for bottom (with argument x/(1-xn))
      if cxn > tol || cxn < -tol
        toGlobal(add, topologyId, mydim, dim-1, cs, df/cxn, x, rf*cxn, X)
      else
        toGlobal(add, topologyId, mydim, dim-1, cs, df, x, T(0), X)
      end
      # apply xn times the tip
      corner = get(cs)
      X .+= (rf*xn) .* corner
      cs.state += 1
    end
  end
end



function jacobianTransposed(add::Bool, topologyId::UInt32, mydim::Integer, dim::Integer, cs::CornerStorage{T}, df::T,
                            x::AbstractVector{T}, rf::T, jt::AbstractArray{T,2}) where {T<:Real}
  rows = Base.size(jt,1)
  cdim = Base.size(jt,2)
  @argcheck rows >= dim
  @argcheck mydim >= dim
  @argcheck length(x) == mydim
  if dim == 0
    cs.state += 1
  else
    xn = df*x[dim]
    cxn = T(1) - xn

    cs2 = CornerStorage(cs.corners, cs.state)

    if TypesImpl.isPrism(topologyId, mydim, mydim-dim)
      # apply (1-xn) times Jacobian for bottom
      jacobianTransposed(add, topologyId, mydim, dim-1, cs2, df, x, rf*cxn, jt)
      # apply xn times Jacobian for top
      jacobianTransposed(true, topologyId, mydim, dim-1, cs2, df, x, rf*xn, jt)
      # compute last row as difference between top value and bottom value
      toGlobal(add, topologyId, mydim, dim-1, cs, df, x, -rf, Base.view(jt,dim,:))
      toGlobal(true, topologyId, mydim, dim-1, cs, df, x, rf, Base.view(jt,dim,:))
    else
      @assert TypesImpl.isPyramid(topologyId, mydim, mydim-dim)

      # In the pyramid case, we need a transformation Tb: B -> R^n for the
      # base B \subset R^{n-1}. The pyramid transformation is then defined as
      #   T: P \subset R^n  -> R^n
      #      (x, xn)        |-> (1-xn) Tb(x*) + xn t  (x \in R^{n-1}, xn \in R)
      # with the tip of the pyramid mapped to t and x* = x/(1-xn)
      # the projection of (x,xn) onto the base.
      #
      # For the Jacobi matrix DT we get
      #   DT = ( A | b )
      # with A = DTb(x*)   (n x n-1 matrix)
      #  and b = dT/dxn    (n-dim column vector).
      # Furthermore
      #   b = -Tb(x*) + t + \sum_i dTb/dx_i(x^*) x_i/(1-xn)
      #
      # Note that both A and b are not defined in the pyramid tip (x=0, xn=1)!
      # Indeed for B the unit square, Tb mapping B to the quadrilateral given
      # by the vertices (0,0,0), (2,0,0), (0,1,0), (1,1,0) and t=(0,0,1), we get
      #
      #   T(x,y,xn) = ( x(2-y/(1-xn)), y, xn )
      #               / 2-y/(1-xn)  -x   0 \
      #  DT(x,y,xn) = |    0         1   0 |
      #               \    0         0   1 /
      # which is not continuous for xn -> 1, choose for example
      #   x=0,    y=1-xn, xn -> 1   --> DT -> diag(1,1,1)
      #   x=1-xn, y=0,    xn -> 1   --> DT -> diag(2,1,1)
      #
      # However, for Tb affine-linear, Tb(y) = My + y0, DTb = M:
      #   A = M
      #   b = -M x* - y0 + t + \sum_i M_i x_i/(1-xn)
      #     = -M x* - y0 + t + M x*
      #     = -y0 + t
      # which is continuous for xn -> 1. Note that this b is also given by
      #   b = -Tb(0) + t + \sum_i dTb/dx_i(0) x_i/1
      # that is replacing x* by 1 and 1-xn by 1 in the formular above.
      #
      # For xn -> 1, we can thus set x*=0, "1-xn"=1 (or anything != 0) and get
      # the right result in case Tb is affine-linear.

      tol = 16 * eps(T)

      # The second case effectively results in x* = 0
      dfcxn = (cxn > tol || cxn < -tol) ? T(df / cxn) : T(0)

      # initialize last row
      # b =  -Tb(x*)
      # (b = -Tb(0) = -y0 in case xn -> 1 and Tb affine-linear)
      toGlobal(add, topologyId, mydim, dim-1, cs, dfcxn, x, -rf, jt[dim,:])
      # b += t
      jt[dim,:] .+= rf .* get(cs)
      cs.state += 1
      # apply Jacobian for bottom (with argument x/(1-xn)) and correct last row
      if add
        jt2 = zeros(T, dim-1, cdim)
        # jt2 = dTb/dx_i(x*)
        jacobianTransposed(false, topologyId, mydim, dim-1, cs2, dfcxn, x, rf, jt2)
        # A = dTb/dx_i(x*)                      (jt[j], j=0..dim-1)
        # b += \sum_i dTb/dx_i(x*) x_i/(1-xn)   (jt[dim-1])
        # (b += 0 in case xn -> 1)
        for j = 1:dim-1
          jt[j,:] .+= jt2[j,:]
          jt[dim,:] .+= dfcxn*x[j] .* jt2[j,:]
        end
      else
        # jt = dTb/dx_i(x*)
        jacobianTransposed(false, topologyId, mydim, dim-1, cs2, dfcxn, x, rf, jt)
        # b += \sum_i dTb/dx_i(x*) x_i/(1-xn)
        for j = 1:dim-1
          jt[dim,:] .+= dfcxn*x[j] .* jt[j,:]
        end
      end
    end
  end
end


function affine(topologyId::UInt32, mydim::Integer, dim::Integer, cs::CornerStorage{T}, jt::Array{T,2})::Bool where {T<:Real}
  if dim == 0
    cs.state += 1
    return true
  else
    orgBottom = get(cs)
    if !affine(topologyId, mydim, dim-1, cs, jt)
      return false
    end

    orgTop = get(cs)
    tol = 16 * eps(T)

    if TypesImpl.isPrism(topologyId, mydim, mydim-dim)
      jtTop = zeros(T, Base.size(jt,1), Base.size(jt,2))
      if !affine(topologyId, mydim, dim-1, cs, jtTop)
        return false
      end

      # check whether both jacobians are identical
      norm = T(0)
      for i = 1:dim-1
        norm += sum((jtTop[i,:] .- jt[i,:]).^2)
      end
      if norm >= tol
        return false
      end
    else
      cs.state += 1
    end
    jt[dim,:] .= orgTop .- orgBottom
    return true
  end
end

end # module MultiLinearGeometryImpl
