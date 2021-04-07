# Gas Dynamic Equations for Normal Shock Waves

# A struct called StreamConditions which stores the state of a particular flow
# Functions for generating the normal shock ratios (functions only of Mach number for a normal shock wave)
# Function for generating new flow condition struct for downstream flow, for given upstream StreamConditions struct

# Thermodynamic parameters
gamma = 1.4

struct StreamConditions
    MachNumber
    Temperature
    Pressure
    Density
end

FreeStream = StreamConditions(2,350,0.6,0.1225)


function M2(M1)

    return (  (M1^2*(gamma-1) + 2)  /  (2*gamma*M1^2 - (gamma-1)) )^0.5

end


function P2(M1)

    return (2*gamma*M1^2 - (gamma -1)) / (gamma+1)

end

function T2(M1)

    return (2*gamma*M1^2 - (gamma-1))*((gamma-1)*M1^2+2) / ((gamma+1)^2*M1^2)

end

function RHO2(M1)

    return ((gamma+1)*M1^2)/((gamma-1)*M1^2+2)

end 



function NormalShock(flow)

    m2   = M2(flow.MachNumber)
    t2   = flow.Temperature*T2(flow.MachNumber)
    p2   = flow.Pressure*P2(flow.MachNumber)
    rho2 = flow.Density*RHO2(flow.MachNumber)
    
    downstream_flow = StreamConditions(m2,t2,p2,rho2)

    return downstream_flow

end


ds = NormalShock(FreeStream)

println(ds.MachNumber)
println(ds.Temperature)
println(ds.Density)
