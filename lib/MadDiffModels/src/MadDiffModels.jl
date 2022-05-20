"""
    MadDiffModels

`MadDiffModels` is a submodule of `MadDiff`. `MadDiffModels` allows modeling nonlinear optimization problem of the following form:
```
minimize:   f(x)
subject to: xl ≤   x  ≤ xu
            gl ≤ g(x) ≤ gu,
```
where:
- `x ∈ R^n` is the decision variable.
- `f : R^n -> R`   is the objective function
- `g : R^n -> R^m` is the constraint mapping.

The model is constructed as an `NLPModel` (see https://github.com/JuliaSmoothOptimizers/NLPModels.jl), and can be solved with various NLP solvers such as:
- MadNLP (https://github.com/MadNLP/MadNLP.jl)
- Ipopt (https://github.com/JuliaSmoothOptimizers/NLPModelsIpopt.jl)
- Knitro (https://github.com/JuliaSmoothOptimizers/NLPModelsKnitro.jl)
"""
module MadDiffModels

import Base: ==, >=, <=, <, >, getproperty, getindex, setindex!, add_sum, copyto!, ndims, size, copyto!, convert, string, print, show, log2, max, csch, acoth, tanh, asecd, ^, *, cospi, csc, abs, hypot, /, rem2pi, -, acos, min, exp10, tan, acsc, log10, log1p, \, acotd, sin, -, exp2, acot, acosh, atan, mod, cotd, atanh, sind, rem, mod2pi, secd, atan, sinh, +, cot, cosh, acsch, sinpi, atand, cos, acosd, rad2deg, +, tand, cscd, asinh, cbrt, asin, transpose, inv, acscd, expm1, log, deg2rad, asec, asind, cosd, abs2, sech, sqrt, asech, sec, exp, coth

import MadDiffCore
import NLPModels: NLPModels, jac_structure!, hess_structure!, obj, cons!, grad!, jac_coord!, hess_coord!, get_nvar, get_ncon, get_nnzh, get_nnzj, get_x0, get_y0, get_lvar, get_uvar, get_lcon, get_ucon

include("model.jl")
include("api.jl")

end # module
