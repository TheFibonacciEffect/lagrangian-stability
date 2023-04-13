using Zygote
using Plots
using Rotations

# define the state of the system
cube_sidelength = 1.0
intertia_tensor_cube = [1/6*cube_sidelength^2 0 0; 0 1/6*cube_sidelength^2 0; 0 0 1/6*cube_sidelength^2]
mass_cube = 1.0
μ_cube = [1.0, 0.0, 0.0]
θ_cube = 0.0

# inital parameters

