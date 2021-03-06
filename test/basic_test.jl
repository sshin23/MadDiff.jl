x = Variable("x")
expr = sin(x[2]^4) + 2

@test func(x[2]) == func(copy(x[2]))
@test deriv(x[2]) == deriv(copy(x[2]))
@test func(expr) == func(copy(expr))
@test deriv(expr) == deriv(copy(expr))
