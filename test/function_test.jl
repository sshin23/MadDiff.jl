for (forig,dorig,points) in [
    (x->one(x[1]), Dict(), [randn(1) for i=1:3]),
    (x->zero(x[1]), Dict(), [randn(1) for i=1:3]),
    (x->+x[1], Dict(1=>x->1.), [randn(1) for i=1:3]),
    (x->x[1]+1., Dict(1=>x->1.), [randn(1) for i=1:3]),
    (x->2x[1], Dict(1=>x->2.), [randn(1) for i=1:3]),
    (x->-x[1], Dict(1=>x->-1.), [randn(1) for i=1:3]),
    (x->x[1]^2, Dict(1=>x->2x[1]), [randn(1) for i=1:3]),
    (x->x[1]^2.5, Dict(1=>x->2.5x[1]^1.5), [rand(1) for i=1:3]),
    (x->1/x[1], Dict(1=>x->-x[1]^-2), [randn(1) for i=1:3]),
    (x->sin(x[1]), Dict(1=>x->cos(x[1])), [randn(1) for i=1:3]),
    (x->erf(x[1]), Dict(1=>x->(2/sqrt(pi)) * exp(-x[1]^2)), [randn(1) for i=1:3]),
    (x->cos(x[1]), Dict(1=>x->-sin(x[1])), [randn(1) for i=1:3]),
    (x->exp(x[1]), Dict(1=>x->exp(x[1])), [randn(1) for i=1:3]),
    (x->x[1]+x[2], Dict(1=>x->1.,2=>x->1.), [randn(2) for i=1:3]),
    (x->x[1]-x[2], Dict(1=>x->1.,2=>x->-1.),[randn(2) for i=1:3]),
    (x->x[1]*x[2], Dict(1=>x->x[2],2=>x->x[1]),[randn(2) for i=1:3]),
    (x->x[1]/x[2], Dict(1=>x->1/x[2],2=>x->-x[1]/x[2]^2),[randn(2) for i=1:3]),
    (x->x[1]^x[2], Dict(1=>x->1/x[2],1=>x->x[2]*x[1]^(x[2]-1),2=>x->x[1]^x[2]*log(x[1])),[rand(2) for i=1:3]),
    (x->sin(x[1]/x[2]), Dict(1=>x->cos(x[1]/x[2])/x[2],2=>x->-cos(x[1]/x[2])*x[1]/x[2]^2), [randn(2) for i=1:3]),
    (x->cos(x[1]-x[2])^2, Dict(1=>x->-2cos(x[1]-x[2])sin(x[1]-x[2]),2=>x->2cos(x[1]-x[2])sin(x[1]-x[2])), [randn(2) for i=1:3]),
    (x->exp(2^x[1]+x[2]^3), Dict(1=>x->exp(2^x[1]+x[2]^3)*2^x[1]*log(2) ,2=>x->exp(2^x[1]+x[2]^3)*3*x[2]^2), [randn(2) for i=1:3]),
    (x->x[1]/log(x[2]^2+9.), Dict(1=>x->1/log(x[2]^2+9.), 2=>x->-2x[2]x[1]/(x[2]^2+9.)/log(x[2]^2+9.)^2), [randn(2) for i=1:3]),
    (x->sum(cos(x[i]) for i=1:4), Dict(i=>x->-sin(x[i]) for i=1:4), [randn(4) for i=1:3]),
    (x->beta(x[1],1), Dict(1=>x->beta(x[1],1)*(digamma(x[1])-digamma(x[1]+1))), [randn(1) for i=1:3]),
    (x->beta(1,x[1]), Dict(1=>x->beta(1,x[1])*(digamma(x[1])-digamma(1+x[1]))), [randn(1) for i=1:3]),
]
    expr = forig(Variable())
    f = func(expr)
    @test compare(forig,f, points)
    
    for (i,d) in deriv(expr)
        @test compare(dorig[i],d,points)
    end
end

for (forig,dorig,xpoints,ppoints) in [
    ((x,p)->x[1]+p[1], Dict(1=>(x,p)->1.), [randn(1) for i=1:3],[randn(1) for i=1:3]),
    ((x,p)->x[1]^p[1], Dict(1=>(x,p)->p[1]*x[1]^(p[1]-1)), [rand(1) for i=1:3],[randn(1) for i=1:3]),
    ((x,p)->x[1]/p[1], Dict(1=>(x,p)->1/p[1]),[randn(1) for i=1:3],[randn(1) for i=1:3]),
    ((x,p)->beta(x[1],p[1]), Dict(1=>(x,p)->beta(x[1],p[1])*(digamma(x[1])-digamma(x[1]+p[1]))), [randn(1) for i=1:3], [randn(1) for i=1:3]),
    ((x,p)->beta(p[1],x[1]), Dict(1=>(x,p)->beta(p[1],x[1])*(digamma(x[1])-digamma(p[1]+x[1]))), [randn(1) for i=1:3], [randn(1) for i=1:3])
]
    expr = forig(Variable(),Parameter())
    f = func(expr)
    @test compare(forig,f, xpoints,ppoints)
    
    for (i,d) in deriv(expr)
        @test compare(dorig[i],d,xpoints,ppoints)
    end
end

for (forig,xpoints,ppoints) in [
    ((x,p)->sum(1 +beta(sin(x[1]+erf(p[2])/p[1])^2,cos(p[1]*x[4])/x[5]^2-x[3])+sin(x[i]) for i=1:7),[randn(7) for i=1:3], [randn(2) for i=1:3]),
    ((x,p)->1 - p[1] - sum(sin(sin(sin(cos(beta(1,x[i])+p[2])-1.)+4)^9)/abs2(x[3]) for i=1:7),[randn(7) for i=1:3], [randn(2) for i=1:3]),
]
    expr = forig(Variable("xx"),Parameter("pp"))
    f = func(expr)
    @test compare(forig,f, xpoints,ppoints)
end
