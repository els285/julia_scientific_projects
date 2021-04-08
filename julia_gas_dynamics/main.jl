include("GasDy.jl")
using .GasDynamics

f1 = GasDynamics.BaseFluid.Fluid("air",1.4,287)

flow = GasDynamics.BaseFluid.MakeFlow(f1,2.5;Density=1,Pressure=1e5)

downstream_normal = GasDynamics.Shocks.NormalShock.ApplyNormalShock(flow)

println(downstream_normal.Pressure/(downstream_normal.Temperature*downstream_normal.Density))

downstream_oblique = GasDynamics.Shocks.ObliqueShock.ApplyObliqueShock(flow,0.3)

println(downstream_oblique.Pressure/(downstream_oblique.Temperature*downstream_oblique.Density))