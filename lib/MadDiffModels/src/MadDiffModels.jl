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

import Base: show, add_sum, getindex, setindex!
import MadDiffCore
import NLPModels: NLPModels, jac_structure!, hess_structure!, obj, cons!, grad!, jac_coord!, hess_coord!

include("model.jl")
include("api.jl")

export MadDiffModel, variable, constraint, objective, parameter, instantiate!, value, setvalue, set_lower_bound, set_upper_bound, lower_bound, upper_bound, set_lower_bound, set_upper_bound, lower_bound, upper_bound, dual

end # module
