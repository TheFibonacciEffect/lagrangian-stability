using Zygote
using Plots
using Rotations

using ModelingToolkit
# define the state of the system
@constants cube_sidelength = 1.0
@constants intertia_tensor_cube = [1/6*cube_sidelength^2 0 0; 0 1/6*cube_sidelength^2 0; 0 0 1/6*cube_sidelength^2]
@constants mass_cube = 1.0 μ_cube = [1.0, 0.0, 0.0] θ_cube = 0.0
@constants μx=1 μy=0 μz=0
@constants Bx = 1 By = 0 Bz = 0
B = [Bx; By; Bz]
μ = [μx; μy; μz]

@variables t α(t) β(t) γ(t) ω_x(t) ω_y(t) ω_z(t) x(t) y(t) z(t) # independent and dependent variables
@variables dαdt(t) dβdt(t) dγdt(t) dω_xdt(t) dω_ydt(t) dω_zdt(t) dxdt(t) dydt(t) dzdt(t) # derivatives of dependent variables
@variables dxdt(t) dydt(t) dzdt(t) # derivatives of dependent variables
D = Differential(t)

# compute omega
R = RotYXZ(α,β,γ)
RtR = D.(transpose(R)) * R
RtR = expand_derivatives(RtR)
ω_x = -RtR[3,2]
ω_y = RtR[3,1]
ω_z = -RtR[2,1]

ω = [ω_x; ω_y; ω_z]

# define lagrangian
r = [x; y; z]
r_dot = D.(r) .|> expand_derivatives
L = 1/2 * mass_cube * (r_dot' * r_dot) + 1/2 * ω' * R * intertia_tensor_cube * transpose(R) * ω + (R*μ)'*B
L_subst = substitute(L, [D(x) => dxdt, D(y) => dydt, D(z) => dzdt, D(α) => dαdt, D(β) => dβdt, D(γ) => dγdt, D(ω_x) => dω_xdt, D(ω_y) => dω_ydt, D(ω_z) => dω_zdt])

# euler lagrange equations
# position
function get_el_equation(x,dxdt)
    Dxdot = Differential(dxdt)
    Dx = Differential(x)
    LHS = D(Dxdot(L_subst))
    RHS = Dx(L_subst)
    return LHS, RHS    
end

get_el_equation(x,dxdt)
get_el_equation(y,dydt)
get_el_equation(z,dzdt)
get_el_equation(α,dαdt)
get_el_equation(β,dβdt)
get_el_equation(γ,dγdt)

expand_derivatives(get_el_equation(x,dxdt)[1]) # this actually works!
expand_derivatives(get_el_equation(α,dαdt)[1],false) # why is this zero?

# expand_derivatives(get_el_equation(α,dαdt)[2],false) # why does this raise an error?
