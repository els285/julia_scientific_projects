# ode 45
function RK45(y_coefs,x_coefs,boundary,step_size)
    x0,xN,y0 = boundary[1],boundary[2],boundary[3]
    
    function F(x,y,x_coefs,y_coefs)
        # Generalised function for the derivative term
        
        x_sigma=0
        for (i,c) in enumerate(x_coefs)
            x_sigma += c*x^i
        end
        
        y_sigma=0
        for (i,c) in enumerate(y_coefs)
            y_sigma += c*y^i
        end

        return -(x_sigma + y_sigma)
    end
        
    y = y0
    for x in x0:step_size:xN
        k1 = step_size*F(x,y,x_coefs,y_coefs)
        k2 = step_size*F(x + step_size/2 , y + k1/2 , x_coefs, y_coefs)
        k3 = step_size*F(x + step_size/2 , y + k2/2 , x_coefs, y_coefs)
        k4 = step_size*F(x + step_size   , y + k3   , x_coefs, y_coefs)
        y = y + k1/6 + k2/3 + k3/3 + k4/6
    end
    print(y)
    
end

RK45((-2),(-1,3),(0,10,0),0.01)