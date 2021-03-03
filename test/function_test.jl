for (forig,dorig,points) in [
    (x->x[1], Dict(1=>x->1.), [randn(1) for i=1:3]),
    (x->x[1]+1., Dict(1=>x->1.), [randn(1) for i=1:3]),
    (x->2x[1], Dict(1=>x->2.), [randn(1) for i=1:3]),
    (x->-x[1], Dict(1=>x->-1.), [randn(1) for i=1:3]),
    (x->x[1]^2, Dict(1=>x->2x[1]), [randn(1) for i=1:3]),
    (x->x[1]^2.5, Dict(1=>x->2.5x[1]^1.5), [rand(1) for i=1:3]),
    (x->1/x[1], Dict(1=>x->-x[1]^-2), [randn(1) for i=1:3]),
    (x->sin(x[1]), Dict(1=>x->cos(x[1])), [randn(1) for i=1:3]),
    (x->cos(x[1]), Dict(1=>x->-sin(x[1])), [randn(1) for i=1:3]),
    (x->exp(x[1]), Dict(1=>x->exp(x[1])), [randn(1) for i=1:3]),
    (x->x[1]+x[2], Dict(1=>x->1.,2=>x->1.), [randn(2) for i=1:3]),
    (x->x[1]-x[2], Dict(1=>x->1.,2=>x->-1.),[randn(2) for i=1:3]),
    (x->x[1]*x[2], Dict(1=>x->x[2],2=>x->x[1]),[randn(2) for i=1:3]),
    (x->x[1]/x[2], Dict(1=>x->1/x[2],2=>x->-x[1]/x[2]^2),[randn(2) for i=1:3]),
    (x->sin(x[1]/x[2]), Dict(1=>x->cos(x[1]/x[2])/x[2],2=>x->-cos(x[1]/x[2])*x[1]/x[2]^2), [randn(2) for i=1:3]),
    (x->cos(x[1]-x[2])^2, Dict(1=>x->-2cos(x[1]-x[2])sin(x[1]-x[2]),2=>x->2cos(x[1]-x[2])sin(x[1]-x[2])), [randn(2) for i=1:3]),
    (x->exp(2^x[1]+x[2]^3), Dict(1=>x->exp(2^x[1]+x[2]^3)*2^(x[1]-1)*log(2) ,2=>x->exp(2^x[1]+x[2]^3)*3*x[2]^2), [randn(2) for i=1:3]),
    (x->x[1]/log(x[2]^2+9.), Dict(1=>x->1/log(x[2]^2+9.), 2=>x->-2x[2]x[1]/(x[2]^2+9.)/log(x[2]^2+9.)^2), [randn(2) for i=1:3]),
]
    expr = forig(Variable())
    f = func(expr)
    @test compare(forig,f, points)
    
    for (i,d) in deriv(expr)
        @test compare(dorig[i],d,points)
    end
end
