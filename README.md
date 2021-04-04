[![build](https://github.com/sshin23/SimpleNL.jl/actions/workflows/test.yml/badge.svg)](https://github.com/sshin23/SimpleNL.jl/actions/workflows/test.yml)
[![codecov](https://codecov.io/gh/sshin23/SimpleNL.jl/branch/main/graph/badge.svg?token=U6NMMW0IT5)](https://codecov.io/gh/sshin23/SimpleNL.jl)
# SimpleNL.jl
 
SimpleNL.jl is a simple algebraic modeling/differentiation package. The package is under development.

## Installation
```julia
pkg> add https://github.com/sshin23/SimpleNL.jl.git
```

## Usage
### Nonlinear Expressions
SimpleNL.jl provides a flexible user-interface for writing nonlinear expressions and evaluating the expressions and functions. For example,
```julia
using SimpleNL

x = Variable()
p = Parameter()
expr = x[1]^2 + exp(x[2]^p[1])/2 + log(x[3]+p[2])
println(expr) # x[1]^2 + exp(x[2]^p[1])/2 + log(x[3] + p[2])

x0 = [0.,0.5,1.5]
p0 = [2,0.5]

f = function_evaluator(expr)
println("f(x0,p0) = $(f(x0,p0))") # f(x0,p0) = 1.3351598889038159

y0 = zeros(3)
d = gradient_evaluator(expr)
d(y0,x0,p0)
println("d_2(x0,p0) = $y0") # d_2(x0,p0) = [0.0, 0.6420127083438707, 0.5]
```

### Nonlinear Programming
SimpleNL.jl provides a simple user-interface for creating nonlinear prgogramming models and allows solving the created models with optimization solvers (Ipopt and MadNLP.jl). The syntax is as follows:
```julia
using SimpleNL, Ipopt

m = SimpleNL.Model(Ipopt.Optimizer; print_level=3) 

x = [variable(m; start=mod(i,2)==1 ? -1.2 : 1.) for i=1:1000]
objective(m, sum(100(x[i-1]^2-x[i])^2+(x[i-1]-1)^2 for i=2:1000))
for i=1:998
    constraint(m, 3x[i+1]^3+2*x[i+2]-5+sin(x[i+1]-x[i+2])sin(x[i+1]+x[i+2])+4x[i+1]-x[i]exp(x[i]-x[i+1])-3 == 0)
end

instantiate!(m) # optional; this will make the model ready to be solved
optimize!(m)
```

## How it works?
SimpleNLModels.jl constructs first and second derivative functions off-line (i.e., prior to calling the optimization solver) by applying operator overloading-based automatic differentiation on _functions_. The exact _derivative functions_ can be obtained as results. A benchmark code is [here](https://github.com/sshin23/SimpleNLModels.jl/blob/main/benchmark/benchmark.jl).

<img src="/benchmark/output/luksanvlcek.png" width="400"/><img src="/benchmark/output/hehnandrea.png" width="400"/>

## Bug reports and support
Please report issues and feature requests via the [Github issue tracker](https://github.com/sshin23/SimpleNLModels.jl/issues). 
