using ModelingToolkit

@variables t x(t)   # independent and dependent variables
@parameters τ       # parameters
@constants h = 1    # constants have an assigned value
D = Differential(t) # define an operator for the differentiation w.r.t. time

LHS(x) = D(x)
# your first ODE, consisting of a single equation, the equality indicated by ~
@named fol = ODESystem([LHS(x) ~ (h - x) / τ])

using DifferentialEquations: solve

prob = ODEProblem(fol, [x => 0.0], (0.0, 10.0), [τ => 3.0])
# parameter `τ` can be assigned a value, but constant `h` cannot
sol = solve(prob)

using Plots
plot(sol)
