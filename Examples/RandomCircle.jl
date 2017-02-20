module RandomCircle

function rand_normal(density::Int64, mean::Float64, stdev::Float64)
  if stdev <= 0.0
      error("standard deviation must be positive")
  end
  const u1 = rand(density)
  const u2 = rand(density)
  const r = sqrt( -2.0*log(u1) )
  const Θ = 2.0*pi*u2
  mean + stdev*sin(Θ).*r
end

function random_circle(center::Array{Float64,1}, radius::Float64, density::Int64)
  radii = rand_normal(density, radius, 0.03)
  # radii = ones(density)
  # angles = rand(density)*2*pi
  angles = linspace(0.,1.-1/density,density)*2*pi
  center_x = center[1]
  center_y = center[2]
  [[radii[i]*cos(angles[i]) + center_x, radii[i]*sin(angles[i]) + center_y] for i in 1:density]
end
end
