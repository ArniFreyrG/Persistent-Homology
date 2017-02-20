using RandomCircle: random_circle
using VietorisRipsComplex: VRComplex
using BasicHomology: compute_homology
using VietorisRipsFiltration: VRFiltration
using PIntervals: p_intervals
using PyPlot, PyCall
@pyimport matplotlib.patches as patch

const stddev = 0.5
const point_density = 15
const circle = random_circle([0.0,0.0], stddev, point_density)

const point_list = circle
const ɛ_min = 0.0
const ɛ_max = 0.5
const filtration_density = 10000

V = VRFiltration(point_list, ɛ_min, ɛ_max, filtration_density)
P = p_intervals(V)

using PyPlot, PyCall
@pyimport matplotlib.patches as patch
const fig = PyPlot.figure()
ax_1 = fig[:add_subplot](211)
ax_1[:set_aspect]("equal")
for p in point_list
  ax_1[:add_artist](patch.Circle((p[1],p[2]), ɛ_max, ec="none", alpha=0.2))
end
ax_1[:set_xlim]([-1,1])
ax_1[:set_ylim]([-1,1])
const X = [point[1] for point in point_list]
const Y = [point[2] for point in point_list]
PyPlot.plot(X, Y, marker=".", linestyle="None")
ax_2 = fig[:add_subplot](212)
i = 30
axhline(y=i, linewidth=3, color=[0,0,0,1])
i -= 1
for dim in P
  for p in dim
    if p[1] != p[2]
      PyPlot.plot(ɛ_max*p/filtration_density,[i,i], "k-")
      i -= 1
    end
  end
  axhline(y=i, linewidth=3, color=[0,0,0,1])
  i -= 1
end
PyPlot.show()
