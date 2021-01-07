#==================
This file defines two function.

# Wave Drag
The first estimates the wave drag for a body of revolution with fixed length and diameter.
The length, and diameter are hardcoded for this problem.
    cd, x, r = wavedrag(rpts)

# Inputs
rpts: a vector of 6 radii at fixed locations between the nose and maximum diameter (5) of the body

# Outputs
cd: the drag coefficient (normalized by the max cross sectional area)
x: a dense sampling of points from 0 to the 100 (the body length) that defines the geometry
r: the corresponding radii at those x locations.

# Example:
rpts = [0.1, 0.198, 0.7, 1.8, 3.200, 4.550]
cd, x, r = wavedrag(rpts)

println("cd = ", cd)

using PyPlot
figure()
plot(x, r)

# Minimum Wave Drag
The second function returns the minimum drag solution for this fixed length and diameter, 
an analytic solution derived by Sears and Haack.  This is provided just for comparison.

    cd, x, r = searshaack()

# Outputs
cd: the minimum drag coefficient for this length and diameter
x: corresponding shape (axial x locations)
r: corresponding shape (radii)

# Example:
cdmin, xmin, rmin = searshaack()

println("cdmin = ", cdmin)

using PyPlot
figure()
plot(xmin, rmin)
===================#


using FLOWMath  # you will need to install from package manager: Pkg.add("FLOWMath")

# ------- add functionality to Akima in FLOWMath for second derivative ---------
function second_deriv(spline::Akima, x)

    j = FLOWMath.findindex(spline.xdata, x)

    # evaluate polynomial
    dx = x - spline.xdata[j]
    d2ydx2 = 2*spline.p2[j] + 6*spline.p3[j]*dx

    return d2ydx2
end

second_deriv(spline::Akima, x::AbstractVector) = second_deriv.(Ref(spline), x)
# -------------------------------------------------------------------------------

# ----------- wave drag function --------------------

function wavedrag(rvar)
    
    # setup geometry with a spline
    L = 100.0
    xpt = L*[0.0, 0.005, 0.01, 0.025, 0.1, 0.2, 0.35, 0.5, 0.65, 0.8, 0.9, 0.975, 0.99, 0.995, 1.0]
    rpt = [0.0; rvar; 5; rvar[end:-1:1]; 0.0]
    spline = Akima(xpt, rpt)

    # setup a finer grid
    x = [0:0.01:0.5; 0.6:98.5; 99.5:0.01:100]
    n = length(x)

    # evaluate second derivative of area
    r = spline(x)
    drdx = gradient(spline, x)
    d2rdx2 = second_deriv(spline, x)
    d2Adx2 = @. 2*pi*(r*d2rdx2 + drdx^2)  # dA^2/dx^2 - second derivative of area

    # setup integrand
    function intg(i, j) 
        if i == j
            return 0.0
        else
            return d2Adx2[i]*d2Adx2[j] * log(abs(x[i] - x[j]))
        end
    end

    # evaluate integral
    I = 0.0
    for i = 1:n-1
        for j = 1:n-1
            
            dx1 = x[i+1] - x[i]
            dx2 = x[j+1] - x[j]
            # trapezoidal rule
            I += (intg(i, j) + intg(i+1, j) + intg(i, j+1) + intg(i+1, j+1))/4 * dx1*dx2 
        end
    end

    # drag coefficient
    Sref = pi*5^2
    cd = -1.0/(2*pi*Sref)*I

    return cd, x, r
end

function searshaack()

    d = 10.0
    L = 100.0

    cd = (pi*d/L)^2

    xSH = .01:.01:1
    rSH = @. d/2*sqrt(sqrt(1-xSH^2)-xSH^2*log((1+sqrt(1-xSH^2))/xSH))
    
    xSH = L/2 .* [-reverse(xSH); xSH] .+ L/2
    rSH = [reverse(rSH); rSH]

    return cd, xSH, rSH
end

# ---------------------------------------------
