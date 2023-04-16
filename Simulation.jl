using Zygote
using Plots
using Rotations

using ModelingToolkit
using DifferentialEquations

# define the state of the system
@constants a = 1.0 #sidelength of cube
I_cube = [1//6*a^2 0 0; 0 1//6*a^2 0; 0 0 1//6*a^2]
@constants mass_cube = 1.0 μ_cube = [1.0, 0.0, 0.0] θ_cube = 0.0
@constants μx=1 μy=0 μz=0
@constants Bx = 1 By = 0 Bz = 0
B = [Bx; 0; 0]
μ = [μx; 0; 0]

@variables t α(t) β(t) γ(t) ω_x(t) ω_y(t) ω_z(t) x(t) y(t) z(t) # independent and dependent variables
@variables dαdt(t) dβdt(t) dγdt(t) dω_xdt(t) dω_ydt(t) dω_zdt(t) dxdt(t) dydt(t) dzdt(t) # derivatives of dependent variables
@variables dxdt(t) dydt(t) dzdt(t) # derivatives of dependent variables
D = Differential(t)

# compute omega
# R = RotYXZ(α,β,γ)
R = RotZ(γ)
RtR = D.(transpose(R)) * R

RtR = expand_derivatives.(RtR,false)
ω_x = -RtR[3,2]
ω_y = RtR[3,1]
ω_z = -RtR[2,1]

ω = [ω_x; ω_y; ω_z]

M = Array{Any}(undef,3,3)
for (i,θ) in enumerate([α, β, γ])
    M[:,i] = @. expand_derivatives(Differential(D(θ))(ω))
end

M*[D(α); D(β); D(γ)]
ω
# define lagrangian
r = [x; y; z]
r_dot = D.(r) .|> expand_derivatives
L = 1/2 * mass_cube * (r_dot' * r_dot) + 1/2 * ω' * I_cube * ω + (R*μ)'*B
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

LHS,RHS = expand_derivatives.(get_el_equation(γ,dγdt))

equations = [LHS ~ RHS,  D(γ)~dγdt]
el_equations = ODESystem(equations, t, [γ, dγdt], [] , name=:el_equations)

# structural_simplify
# this has three variables:
# States (3):
# γ(t) = angle around z axis
# dγdt(t) = angular velocity around z axis
# γˍtt(t) = angular acceleration around z axis
el_simplefied = structural_simplify(el_equations)

# solve the equations
prob = ODEProblem(el_simplefied, [π/2, 0.0,0.0], (0.0, 10.0))

sol = solve(prob)

# plot the solution
# plot(sol, vars=(t, γ, dγdt))
plot(sol, vars=(0,1))

# compare with solution

sysalg2 = ODESystem([D(γ)~dγdt,D(dγdt)~-Bx*μx/I_cube[1,1]*sin(γ)], t, [γ, dγdt], [] , name=:el_equations)
prob2 = ODEProblem(sysalg2, [π/2, 0.0], (0.0, 10.0))
sol2 = solve(prob2)
plot!(sol2, vars=(0,1))

# they match!
