module BasicHomology

abstract SimplicalComplex

function compute_faces(σ::Array{Int64,1})
  k = length(σ)
  return map(i -> [σ[1:(i-1)] ; σ[(i+1):k]], k:-1:1)::Array{Array{Int64,1}}
end

function boundary_matrices(K::SimplicalComplex)
  const simplices = K.simplices
  const D = length(simplices)-1
  matrices = Array{Array{Int64,2},1}(D)
  for i in 1:D
    n = length(simplices[i+1])
    m = length(simplices[i])
    δ_i = zeros(Int64,n,m)
    if n != 0 && m != 0
      for (j,σ) in enumerate(simplices[i+1])
        faces = compute_faces(σ)
        face_indices = map(ρ -> searchsortedfirst(simplices[i], ρ, lt=lexless), faces)
        for (k,face_index) in enumerate(face_indices)
          δ_i[j, face_index] = (-1)^(k%2)
        end
      end
    end
    matrices[i] = δ_i
  end
  matrices
end

function compute_homology(K::SimplicalComplex)
  println("Computing homology")
  tic()
  const D = length(K.simplices)
  const matrices = boundary_matrices(K)
  const ranks = [[0] ; map(δ -> rank(δ), matrices) ; [0]]
  const added_ranks = map(i -> ranks[i] + ranks[i+1], 1:D)
  const lengths = map(S -> length(S), K.simplices[1:length(K.simplices)])
  const H = lengths - added_ranks
  toc()
  for (d,b) in enumerate(H)
    println("Betti number $(d-1) is $b.")
  end
end

end
