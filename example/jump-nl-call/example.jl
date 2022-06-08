# Motivation:
# - need for function-based nonlinear expressions
# - performance (fast/non-allocating callbacks)
# - sparse jacobian/hessian (w/o coloring)
# - symbolic diff can be an alternative, but want to AD if possible

using MadDiff # loads MadDiffCore, MadDiffSpecialFunctions, MadDiffModels, MadDiffMOI

x = Variable() # "sources" of variable

typeof(x)
typeof(x[2])
typeof(Variable(2))

expr = cos(sin(x[1]*x[3]) + 2) # builds expression tree

# ---------------------------------------------------
# evaluators
# ---------------------------------------------------

f = Evaluator(expr)
g = GradientEvaluator(expr)
h = HessianEvaluator(expr)

x0 = randn(3)
y0 = randn(3)
z0 = randn(3,3)

# non-allocating

@time f(x0) 
@time g(y0,x0); y0
@time h(z0,x0); z0

# ---------------------------------------------------
# sparse evaluators
# ---------------------------------------------------

y1 = randn(2)
z1 = randn(3)

sg = SparseGradientEvaluator(expr)
sh = SparseHessianEvaluator(expr)

sg.sparsity
sh.sparsity

# non-allocating, too

@time sg(y1,x0); y1
@time sh(z1,x0); z1

# ---------------------------------------------------
# NLP example - JuMP-like syntax but function-based
# ---------------------------------------------------
N = 1000

m = MadDiffModel()
x = [variable(m; start = mod(i,2)==1 ? -1.2 : 1.) for i=1:N];
objective(m, sum(100(x[i-1]^2-x[i])^2+(x[i-1]-1)^2 for i=2:N));
for i=1:N-2
    constraint(m, 3x[i+1]^3+2*x[i+2]-5+sin(x[i+1]-x[i+2])sin(x[i+1]+x[i+2])+4x[i+1]-x[i]exp(x[i]-x[i+1])-3 == 0);
end
instantiate!(m)

using NLPModelsIpopt
ipopt(m)

# ---------------------------------------------------
# As a JuMP's AD backend - currently works with MOI#master and JuMP#od/moi-nonlinear
# ---------------------------------------------------

using JuMP, Ipopt

m = JuMP.Model(Ipopt.Optimizer)

@variable(m, x[i=1:N], start=mod(i,2)==1 ? -1.2 : 1.)
@NLobjective(m, Min, sum(100(x[i-1]^2-x[i])^2+(x[i-1]-1)^2 for i=2:N))
@NLconstraint(m, [i=1:N-2], 3x[i+1]^3+2*x[i+2]-5+sin(x[i+1]-x[i+2])sin(x[i+1]+x[i+2])+4x[i+1]-x[i]exp(x[i]-x[i+1])-3 == 0)

optimize!(m; differentiation_backend = MadDiffAD()) # set AD backend

# ---------------------------------------------------
# How it works
# ---------------------------------------------------
# - more similar to symbolic diff in that it creates reusable "evaluator"
# - expression tree or tape doesn't need to be recreated every time -- evaluation is non-allocating.
# - key difference from symbolic diff is that we store intermediate calculation results, so more efficient when we have along nested expressions

function sin100(x)
    for i=1:100
        x = sin(x)
    end
    x
end

x = MadDiff.Variable()
expr = sin100(x[1])


x0 = [1.]
y0 = [0.]
g = GradientEvaluator(expr)

# if we do symbolic diff
using Symbolics

Symbolics.@variables x
expr = sin100(x)

d = expand_derivatives(Differential(x)(expr))

using ReverseDiff

ff = x->sin100(x[1])

@time g(y0,x0) # fast
@time ReverseDiff.gradient!(y0,ff,x0); # fast
@time substitute(d,Dict(x=>1.)); # slow

# ---------------------------------------------------
# how summation/"field" works
# ---------------------------------------------------
# the goal is minimizing the type inference

x = MadDiff.Variable()

expr = sum(
    if mod(i,3) == 0
        sin(x[i])
    elseif mod(i,3) == 1
        cos(x[i])
    elseif mod(i,3) == 2
        exp(x[i])
    end
    for i=1:100)

v = Field()

for i=1:100
    if mod(i,3) == 0
        v[i] = sin(x[i])
    elseif mod(i,3) == 1
        v[i] = cos(x[i])
    elseif mod(i,3) == 2
        v[i] = exp(x[i])
    end
end

# ---------------------------------------------------
# Known limitations
# ---------------------------------------------------
# ---------------------------------------------------
# * control flow - cannot handle if
# ---------------------------------------------------

x = Variable(1)
function my_f(x)
    if x >= 0.
        return x^3
    else
        return -x^3
    end
end

# my_f(x) throws an error. so, do instead

function my_f2(x)
    MadDiffCore.ifelse(x >= 0., x^3, -x^3)
end
my_f2(x)

# or register it!
function my_f(x)
    if x >= 0.
        return x^3
    else
        return -x^3
    end
end

function my_d(x)
    if x >= 0.
        return 3*x^2
    else
        return -3*x^2
    end
end

function my_h(x)
    if x >= 0.
        return 6*x
    else
        return -6*x
    end
end

@register_univariate(my_f,my_d,my_h)

expr = my_f(x)
g = GradientEvaluator(expr)
x0 = [1.]
y0 = [0.]
g(y0,x0)

# ---------------------------------------------------
# * heterogeniety in the algebraic structure
# ---------------------------------------------------

x = Variable(1)
randf(x) = rand(MadDiffCore._UNIVARIATE_FUNCTIONS)(x)

expr1 = sum(sin(x) * sin(x) * sin(x) * sin(x) * sin(x) for i=1:100)
expr2 = sum(randf(x) * randf(x) * randf(x) * randf(x) * randf(x) for i=1:100)

x0 = [1.]
y0 = [0.]
g1 = GradientEvaluator(expr1);
g2 = GradientEvaluator(expr2);
@time g1(y0,x0)
@time g2(y0,x0)

# ---------------------------------------------------
# * cannot efficiently handle common expressions
# ---------------------------------------------------
# - when the same expressions appear multiple time, it is better to evaluate the branch once with the accumulated adjoint
# - but MadDiff currently cannot do that
# - not a big issue for classical applications, but

# some ff neural net
W = randn(10,10)
nn_layer(W,x) = tanh.([sum(W[i,j]*x[j] for j=1:size(W,1)) for i=1:size(W,1)])
function nn_network(n,W,x)
    y = x
    for i=1:n
        y = nn_layer(W,y)
    end
    sum(y)
end

x = [MadDiff.Variable(i) for i=1:10]
expr = nn_network(4,W,x);

g = GradientEvaluator(expr)

x0 = randn(10)
y0 = randn(10)

nn4(x) = nn_network(4,W,x)

@time g(y0,x0) # slow
@time ReverseDiff.gradient!(y0,nn4,x0); # fast

# ---------------------------------------------------
# Development roadmap (w/ no specific timeline)
# - more efficient evaluation of common expressions
# - threaded/gpu evaluator
# ---------------------------------------------------

