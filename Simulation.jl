using Zygote
using Plots
using Rotations

using ModelingToolkit
# define the state of the system
@constants cube_sidelength = 1.0
@constants intertia_tensor_cube = [1/6*cube_sidelength^2 0 0; 0 1/6*cube_sidelength^2 0; 0 0 1/6*cube_sidelength^2]
@constants mass_cube = 1.0 μ_cube = [1.0, 0.0, 0.0] θ_cube = 0.0

@variables t α(t) β(t) γ(t) ω_x(t) ω_y(t) ω_z(t) x(t) y(t) z(t) # independent and dependent variables
D = Differential(t)

# compute omega
RtR = D.(transpose(RotYXZ(α,β,γ))) * RotYXZ(α,β,γ)
RtR = expand_derivatives(RtR)
ω_x = -RtR[3,2]
ω_y = RtR[3,1]
ω_z = -RtR[2,1]
using Latexify
latexify(RtR)
latexify(simplify(ω_x))

ω = [ω_x; ω_y; ω_z]

