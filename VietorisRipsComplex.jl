module VietorisRipsComplex
using BasicHomology: SimplicalComplex

# A point is a list of floats
# A simplex is a list of indices (of points)
# A complex is a list of lists of simplices of the same dimension
function compute_simplices(points::Array{Array{Float64,1},1}, simplices::Array{Array{Int64,1},1}, ɛ::Float64)
  next_dim = Array{Array{Int64,1},1}()
  for (i,p) in enumerate(points)
    const filtered = filter(simplex ->
        !any(s -> s==i, simplex) &&                     # We don't want repeated points in simplices
        !any(s -> norm(points[s]-p) >= 2*ɛ, simplex) && # Check all distances are within eps
        !any(s -> s<i, simplex),                        # Don't want to repeat simplices
      simplices)
    if length(filtered) == 0
      continue
    end
    new_simplices = map(simplex -> sort!([simplex ; [i]]), filtered)::Array{Array{Int64,1},1}
    next_dim = sort!([new_simplices ; next_dim], lt=lexless)::Array{Array{Int64,1},1}
  end
  if length(next_dim) > 1
    return [[next_dim] ; compute_simplices(points, next_dim, ɛ)]::Array{Array{Array{Int64,1},1},1}
  else
    return Array{Array{Array{Int64,1},1},1}()
  end
end

type VRComplex <: SimplicalComplex
  simplices::Array{Array{Array{Int64,1},1},1}
  function VRComplex(points::Array{Array{Float64,1}}, ɛ::Float64)
    println("Complex construction")
    tic()
    const zero_dim_simplices = map(i -> [i], 1:length(points))
    simplices = [[zero_dim_simplices] ; compute_simplices(points, zero_dim_simplices, ɛ)]::Array{Array{Array{Int64,1},1},1}
    toc()
    for (d, σ) in enumerate(simplices)
      println(d-1, "-dim: ", length(σ), " simplices")
    end
    new(simplices)
  end
end
end
