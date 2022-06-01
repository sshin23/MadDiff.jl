# # Getting Started
# Automatic 
MadDiff.jl provides a flexible user-interface for writing nonlinear expressions and evaluating the expressions and functions. For example,
```julia
using MadDiff

x = Variable()
p = Parameter()
expr = x[1]^2 + exp(x[2]^p[1])/2 + log(x[3]+p[2])
println(expr) # x[1]^2 + exp(x[2]^p[1])/2 + log(x[3] + p[2])

x0 = [0.,0.5,1.5]
p0 = [2,0.5]

f = function_evaluator(expr)
println("f(x0,p0) = $(f(x0,p0))") # f(x0,p0) = 1.3351598889038159

y0 = zeros(3)
g = gradient_evaluator(expr)
g(y0,x0,p0)
println("g(x0,p0) = $y0") # g(x0,p0) = [0.0, 0.6420127083438707, 0.5]
```

### Nonlinear Programming
MadDiff.jl provides a simple user-interface for creating nonlinear prgogramming models and allows solving the created models using the solvers with `NLPModels.jl` interface (such as [NLPModelsIpopt.jl](https://github.com/JuliaSmoothOptimizers/NLPModelsIpopt.jl) and [MadNLP.jl](https://github.com/MadNLP/MadNLP.jl)). The syntax is as follows:
```julia
using MadDiff, NLPModelsIpopt

m = MadDiffModel(; print_level=3) 

x = [variable(m; start=mod(i,2)==1 ? -1.2 : 1.) for i=1:1000]
objective(m, sum(100(x[i-1]^2-x[i])^2+(x[i-1]-1)^2 for i=2:1000))
for i=1:998
    constraint(m, 3x[i+1]^3+2*x[i+2]-5+sin(x[i+1]-x[i+2])sin(x[i+1]+x[i+2])+4x[i+1]-x[i]exp(x[i]-x[i+1])-3 == 0)
end

instantiate!(m) # this makes the model ready to be solved
ipopt(m)
```

### Use with JuMP
MadDiff.jl can be used as an automatic differentiation backend. The syntax is as follows:
```julia
using MadDiff, JuMP, Ipopt

m = JuMP.Model(Ipopt.Optimizer) 

@variable(m, x[i=1:1000], start=mod(i,2)==1 ? -1.2 : 1.)
@NLobjective(m, Min, sum(100(x[i-1]^2-x[i])^2+(x[i-1]-1)^2 for i=2:1000))
@NLconstraint(m, [i=1:998], 3x[i+1]^3+2*x[i+2]-5+sin(x[i+1]-x[i+2])sin(x[i+1]+x[i+2])+4x[i+1]-x[i]exp(x[i]-x[i+1])-3 == 0)

optimize!(m; differentiation_backend = MadDiffAD())
```
