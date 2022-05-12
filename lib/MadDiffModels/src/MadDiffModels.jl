module MadDiffModels

import Base: ==, >=, <=, getindex, setindex!, add_sum, copyto!, ndims, size, copyto!, convert, string, print, show, log2, max, csch, acoth, tanh, asecd, ^, *, cospi, csc, abs, hypot, /, rem2pi, -, acos, min, exp10, tan, acsc, log10, log1p, \, acotd, sin, -, exp2, acot, acosh, atan, mod, cotd, atanh, sind, rem, mod2pi, secd, atan, sinh, +, cot, cosh, acsch, sinpi, atand, cos, acosd, rad2deg, +, tand, cscd, asinh, cbrt, asin, transpose, inv, acscd, expm1, log, deg2rad, asec, asind, cosd, abs2, sech, sqrt, asech, sec, exp, coth

import MadDiffCore: MadDiffCore, Sink, Field, Gradient, Jacobian, Hessian, Expression, Variable, Parameter, Constant, f_nargs_1, f_nargs_2, ExpressionSum, inner, SparseJacobian, SparseLagrangianHessian, non_caching_eval, Dummy, DUMMY, index

import NLPModels: finalize, AbstractNLPModel, obj, grad!, cons!, jac_coord!, hess_coord!, hess_structure!, jac_structure!, NLPModelMeta, get_nvar, get_ncon, get_minimize, get_x0, get_y0, get_nnzj, get_nnzh, get_lvar, get_uvar, get_lcon, get_ucon, Counters as NLPModelsCounters, increment! # get_zl,get_zu

export MadDiffModel, variable, parameter, objective, constraint, num_variables, num_constraints, value, dual, objective_value, setvalue, set_lower_bound, set_upper_bound, lower_bound, upper_bound, instantiate!

include("model.jl")

end # module
