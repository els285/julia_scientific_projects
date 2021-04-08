module Shocks

"""
The Shocks module contains two sub-modules: NormalShock and ObliqueShock
"""

module NormalShock

    """
    Normal Shock Module - Gas Dynamics

    Ethan Simpson - 8th April 2021

    Selection of functions for generating thermodynamic and MachNumber parameters downstream of some normal shock,
        for some specified upwind flow condition
    """

    include("./freestream.jl")
    using .BaseFluid

    export ApplyNormalShock
        
    function M2(M1,gamma)

        return (  (M1^2*(gamma-1) + 2)  /  (2*gamma*M1^2 - (gamma-1)) )^0.5

    end

    function P2(M1,gamma)

        return (2*gamma*M1^2 - (gamma -1)) / (gamma+1)

    end

    function T2(M1,gamma)

        return (2*gamma*M1^2 - (gamma-1))*((gamma-1)*M1^2+2) / ((gamma+1)^2*M1^2)

    end

    function Rho2(M1,gamma)

        return ((gamma+1)*M1^2)/((gamma-1)*M1^2+2)

    end 



    function ApplyNormalShock(flow)

        m1 = flow.MachNumber
        gamma = flow.FluidObj.Ratio_of_Specific_Heats

        m2   = M2(m1,gamma)
        t2   = flow.Temperature*T2(m1,gamma)
        p2   = flow.Pressure*P2(m1,gamma)
        rho2 = flow.Density*Rho2(m1,gamma)

        # Cannot use the old type as it is of a different Namespace?
        # Therefore have to rebuild the Fluid object here
        # Then generate a new Flow
        ds_fluid = BaseFluid.Fluid(flow.FluidObj.Name,flow.FluidObj.Ratio_of_Specific_Heats,flow.FluidObj.Gas_Constant)
        downstream_flow = BaseFluid.Flow(ds_fluid,m2,t2,p2,rho2)

        return downstream_flow

    end

end

##############################################################################################################

############################################################################################################



module ObliqueShock

    """
    Oblique Shock Module - Gas Dynamics

    Ethan Simpson - 8th April 2021

    Selection of functions for generating thermodynamic and MachNumber parameters downstream of some oblique shock,
        for some specified upwind flow condition
    """

    include("./freestream.jl")
    using .BaseFluid

    export ApplyObliqueShock

    function ConeAngle(M1,beta,gamma)

        # Gives cone angle for input Mach number and shock angle
        tan_theta = 2*cot(beta)*(M1^2*sin(beta)^2-1)/(M1^2*(gamma + cos(2*beta))+2)

        return atan(tan_theta)

    end


    function M2(M1,beta,gamma)

        theta = ConeAngle(M1,beta,gamma)

        return 1/(sin(beta - theta))*((1 + 0.5*(gamma-1)*M1^2*sin(beta)^2)/(gamma*M1^2*sin(beta)^2) - 0.5*(gamma-1))^0.5

    end


    function P2(M1,beta,gamma)

        return 1 + 2*gamma/(gamma+1) *(M1^2*sin(beta)^2 - 1)

    end

    function T2(M1,beta,gamma)

        return (2*gamma*M1^2*sin(beta)^2 - (gamma-1)) * ((gamma-1)*M1^2*sin(beta)^2 +2) / ((gamma+1)^2*M1^2*sin(beta)^2)

    end

    function Rho2(M1,beta,gamma)

        return ((gamma+1)*M1^2*sin(beta)^2)  / ((gamma-1)*M1^2*sin(beta)^2 + 2 )

    end 


    function ApplyObliqueShock(flow,beta)

        m1 = flow.MachNumber
        gamma = flow.FluidObj.Ratio_of_Specific_Heats

        m2   = M2(m1,beta,gamma)
        t2   = flow.Temperature*T2(m1,beta,gamma)
        p2   = flow.Pressure*P2(m1,beta,gamma)
        rho2 = flow.Density*Rho2(m1,beta,gamma)
        
        ds_fluid = BaseFluid.Fluid(flow.FluidObj.Name,flow.FluidObj.Ratio_of_Specific_Heats,flow.FluidObj.Gas_Constant)
        downstream_flow = BaseFluid.Flow(ds_fluid,m2,t2,p2,rho2)

        return downstream_flow

    end

end

end