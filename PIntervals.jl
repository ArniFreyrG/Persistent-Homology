module PIntervals
using VietorisRipsFiltration: VRFiltration

type SimplexItem
  is_marked::Bool
  d::SparseVector{Float64}
end

function p_intervals(K::VRFiltration)
  const m = length(K.simplices)
  # TODO: Initialize only when used
  T = map(i -> SimplexItem(false, sparse([])), 1:m)
  # TODO: Multisets
  L = map(i -> Set{Array{Int,1}}(), 1:K.max_dim)
  for (j, σ_j) in enumerate(K.simplices)
    (deg, dim, d) = σ_j
    # Remove unmarked terms
    for k in d.nzind
      if !T[k].is_marked
        d[k] = 0
      end
    end
    dropzeros!(d)
    # Gaussian elimination
    while nnz(d) > 0
      i = maximum(d.nzind)
      τ = T[i].d

      if nnz(τ) == 0
        break
      end
      d -= (d[i]/τ[i])*τ
    end
    dropzeros!(d)
    if nnz(d) == 0
      T[j].is_marked = true
    else
      i = maximum(d.nzind)
      T[i].d = d
      push!(L[dim-1], [K.simplices[i][1], deg])
    end
  end

  for j in 1:m
    k = K.simplices[j][2]
    if T[j].is_marked && nnz(T[j].d) == 0
      push!(L[k], [K.simplices[j][1], 10000]) # Should be infinity, but nevermind
    end
  end
  for (d,l) in enumerate(L)
    println("Dimension $(d-1): $(length(l)) intervals")
  end
  return L
end
end
