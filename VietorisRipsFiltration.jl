module VietorisRipsFiltration
using VietorisRipsComplex: VRComplex
using BasicHomology: compute_faces
# using Nemo

type VRFiltration
  simplices::Array{Tuple{Int64,Int64,SparseVector{Int64,Int64}}}
  max_dim::Int64
  function VRFiltration(points::Array{Array{Float64,1}}, ɛ_min::Float64, ɛ_max::Float64, density::Int64)
    println("VR-Complex.")
    tic()
    K_max = VRComplex(points, ɛ_max)
    toc()
    const dim = length(K_max.simplices)
    const dimension_indices = cumsum(map(d -> length(d), K_max.simplices), 1)
    const m = dimension_indices[dim]
    const ɛ_grid = linspace(2*ɛ_min, 2*ɛ_max, density)

    println("Filtration.")
    tic()
    # Filtration
    simplices = Array{Tuple{Int64,Int64,SparseVector{Int64,Int64}},1}(m)

    # Speed things up here by constructing both
    # degrees and boundary mappings at same time
    for d in 1:(dim)
      val = map(j -> (-1)^((j-1)%2), 1:(d))
      if d == 1
        for (j,σ) in enumerate(K_max.simplices[d])
          simplices[j] = (0, 1, sparsevec([]))
        end
      elseif d == 2
        M = dimension_indices[d-1]
        degrees = map(σ ->
          searchsortedfirst(ɛ_grid, norm(points[σ[1]]-points[σ[2]])), K_max.simplices[d]
        )::Array{Int64,1}
        for (j,σ) in enumerate(K_max.simplices[d])
          faces = compute_faces(σ)
          face_indices = map(ρ -> searchsortedfirst(K_max.simplices[d-1], ρ, lt=lexless), faces)
          # Degree of 1-simplices is determined by pairwise distances.
          degree = searchsortedfirst(ɛ_grid, norm(points[σ[1]]-points[σ[2]]))
          # For computational purposes simplices are now of the following form:
          #     (degree, dimension, δ-image)
          simplices[M+j] = (degree, d, sparsevec(face_indices, val, m))
        end
      else
        M = dimension_indices[d-1]
        M_prev = dimension_indices[d-2]
        for (j,σ) in enumerate(K_max.simplices[d])
          # We iterate through simplices of dimension d, and
          # find all faces and their indices for each simplex.
          faces = compute_faces(σ)
          face_indices = map(ρ -> M_prev + searchsortedfirst(K_max.simplices[d-1], ρ, lt=lexless), faces)
          # Degree of σ is the highest degree among its faces.
          degree = reduce(max, (map(i -> simplices[i][1], face_indices)))
          simplices[M+j] = (degree, d, sparsevec(face_indices, val, m))
        end
      end
    end
    toc()

    # Sort simplices by degrees
    println("Sorting by degree.")
    tic()
    α = sortperm(simplices, by=i->i[1])
    α_inverse = sortperm(α)
    # Sort simplices by degree
    simplices = simplices[α]
    # Adjust δ-images for permutation
    # This can probably be sped up
    for (j,σ) in enumerate(simplices)
      im = σ[3]
      permuted_image = map(i -> (α_inverse[i]::Int64, im[i]::Int64), im.nzind)::Array{Tuple{Int64,Int64},1}
      indices = map(v -> v[1], permuted_image)::Array{Int64,1}
      vals = map(v -> v[2], permuted_image)::Array{Int64,1}
      simplices[j] = (σ[1], σ[2], sparsevec(indices, vals, m))
    end
    toc()

    new(simplices, dim)
  end
end
end
