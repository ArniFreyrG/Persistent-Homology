using RandomCircle: random_circle
using VietorisRipsComplex: VRComplex
using BasicHomology: compute_homology
using VietorisRipsFiltration: VRFiltration
using PIntervals: p_intervals
using PyPlot, PyCall
@pyimport matplotlib.patches as patch

const stddev = 0.5
const point_density = 15
const filtration_density = 10000
const circle1 = random_circle([0.0,0.0], stddev, point_density)
const circle2 = random_circle([1.0,0.0], stddev, point_density)

const ɛ = 0.25
const point_list = [circle1 ; circle2]

const V = VRComplex(point_list, ɛ)
compute_homology(V)
const fig = PyPlot.figure()
ax_1 = fig[:add_subplot](111)
ax_1[:set_aspect]("equal")
for p in point_list
  ax_1[:add_artist](patch.Circle((p[1],p[2]), ɛ, ec="none", alpha=0.2))
end
ax_1[:set_xlim]([-1,2])
ax_1[:set_ylim]([-1,1])
const X = [point[1] for point in point_list]
const Y = [point[2] for point in point_list]
PyPlot.plot(X, Y, marker=".", linestyle="None")
PyPlot.show()
