# # Getting Started
# ## Automatic Differentiation
# `MadDiff` provides a flexible user-interface for evaluating first/second-order derivatives of nonlinear expressions. In the following example, using `MadDiff`, we will create a function, gradient, and Hessian evaluator of the following function:
# ```math
# f(x) = x_1^2 + e^{(x_2^{p_1})/2} + \log(x_2x_3+p_2),
# ```
# where ``x`` is the variable vector, and ``p`` is the parameter vector.

# We first import `MadDiff`.

using MadDiff

# First, we create a `Source` of `Variable`'s.

x = Variable()

# The `Base.getindex!` function is extended so that `x[i]` for any `i` creates an expression for ``x_i``. For example,

x[2]

# We can do a similar thing for `Parameter`'s.

p = Parameter()
p[1]

# Now, we create the nonlienar expression expression.

expr = x[1]^2 + exp(x[2]^p[1])/2 + log(x[2]*x[3]+p[2])

# The function evaluator of the above expression can be created by using `MadDiff.function_evaluator` as follows:

f = function_evaluator(expr)

# Now for a given variable and parameter values, the function can be evaluated as follows.
x0 = [0.,0.5,1.5]
p0 = [2,0.5]
f(x0,p0)

# The gradient evaluator can be created as follows:

y0 = similar(x0)
g = gradient_evaluator(expr)
g(y0,x0,p0)
y0

# The Hessian evaluator can be created as follows:
z0 = zeros(3,3)
h = hessian_evaluator(expr)
h(z0,x0,p0)
z0
# Note that only lower-triangular entries are evaluated.

# The evaluator can be constructed in a sparse format:
sh,ij = sparse_hessian_evaluator(expr);
z1 = zeros(length(ij))
sh(z1,x0,p0)
z1

# The sparse coordinates are:
ij


# ## Nonlinear Programming
# ### Built-in API
# MadDiff provides a built-in API for creating nonlinear prgogramming models and allows solving the created models using NLP solvers (in particular, those that are interfaced with `NLPModels`, such as [NLPModelsIpopt](https://github.com/JuliaSmoothOptimizers/NLPModelsIpopt.jl) and [MadNLP](https://github.com/MadNLP/MadNLP.jl)). We now use `MadDiff`'s bulit-in API to model the following nonlinear program:
# ```math
# \begin{aligned}
# \min_{\{x_i\}_{i=0}^N} &\sum_{i=2}^N  100(x_{i-1}^2-x_i)^2+(x_{i-1}-1)^2
# \text{s.t.} &  3x_{i+1}^3+2*x_{i+2}-5+sin(x_{i+1}-x_{i+2})sin(x_{i+1}+x_{i+2})+4x_{i+1}-x_i e^{x_i-x_{i+1}}-3 = 0
# \end{aligned}
# ```
# We model the problem with:
N = 10000

# First, we create a `MadDiffModel`.
m = MadDiffModel() 

# The variables can be created as follows:
x = [variable(m; start=mod(i,2)==1 ? -1.2 : 1.) for i=1:N];

# The objective can be set as follows:
objective(m, sum(100(x[i-1]^2-x[i])^2+(x[i-1]-1)^2 for i=2:N));

# The constraints can be set as follows:
for i=1:N-2
    constraint(m, 3x[i+1]^3+2*x[i+2]-5+sin(x[i+1]-x[i+2])sin(x[i+1]+x[i+2])+4x[i+1]-x[i]exp(x[i]-x[i+1])-3 == 0);
end

# The important last step is instantiating the model. This step must be taken before calling optimizers.
instantiate!(m)

# To solve the problem with `Ipopt`,
using NLPModelsIpopt
ipopt(m);


# ### MadDiff as a AD backend of JuMP
# MadDiff can be used as an automatic differentiation backend of JuMP. The problem above can be modeled in `JuMP` and solved with `Ipopt` along with `MadDiff`

using JuMP, Ipopt

m = JuMP.Model(Ipopt.Optimizer) 

@variable(m, x[i=1:N], start=mod(i,2)==1 ? -1.2 : 1.)
@NLobjective(m, Min, sum(100(x[i-1]^2-x[i])^2+(x[i-1]-1)^2 for i=2:N))
@NLconstraint(m, [i=1:N-2], 3x[i+1]^3+2*x[i+2]-5+sin(x[i+1]-x[i+2])sin(x[i+1]+x[i+2])+4x[i+1]-x[i]exp(x[i]-x[i+1])-3 == 0)

optimize!(m; differentiation_backend = MadDiffAD())
