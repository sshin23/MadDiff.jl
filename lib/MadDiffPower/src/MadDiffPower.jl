module MadDiffPower

import MadDiffCore

# const PF1 =
MadDiffCore.Expression2{typeof(-), Variable, MadDiffCore.Expression2{typeof(+), MadDiffCore.Expression2{typeof(+), MadDiffCore.Expression2{typeof(*), Float64, MadDiffCore.Expression2{typeof(^), Variable, Float64}}, MadDiffCore.Expression2{typeof(*), Float64, MadDiffCore.Expression2{typeof(*), MadDiffCore.Expression2{typeof(*), Variable, Variable}, MadDiffCore.Expression1{typeof(cos), MadDiffCore.Expression2{typeof(-), Variable, Variable}}}}}, MadDiffCore.Expression2{typeof(*), Float64, MadDiffCore.Expression2{typeof(*), MadDiffCore.Expression2{typeof(*), Variable, Variable}, MadDiffCore.Expression1{typeof(sin), MadDiffCore.Expression2{typeof(-), Variable, Variable}}}}}}

MadDiffCore.Expression2{typeof(-), Variable, MadDiffCore.Expression2{typeof(+), MadDiffCore.Expression2{typeof(+), MadDiffCore.Expression2{typeof(*), Float64, MadDiffCore.Expression2{typeof(^), Variable, Float64}}, MadDiffCore.Expression2{typeof(*), MadDiffCore.Expression2{typeof(*), MadDiffCore.Expression2{typeof(*), Float64, Variable}, Variable}, MadDiffCore.Expression1{typeof(cos), MadDiffCore.Expression2{typeof(-), Variable, Variable}}}}, MadDiffCore.Expression2{typeof(*), MadDiffCore.Expression2{typeof(*), MadDiffCore.Expression2{typeof(*), Float64, Variable}, Variable}, MadDiffCore.Expression1{typeof(sin), MadDiffCore.Expression2{typeof(-), Variable, Variable}}}}}

const F64 = 1.0
x[1] - (F64 *x[2]^2.0 + F64*(x[2]*x[3])*cos(x[4] - x[5]) + 18.861786098824915*(x[2]*x[3])*sin(x[4] - x[5]))



end # module
