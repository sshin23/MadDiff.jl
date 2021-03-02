for (f,points) in [
    (x->x[1], [randn(1) for i=1:3]),
    (x->x[1]+1., [randn(1) for i=1:3]),
    (x->2x[1], [randn(1) for i=1:3]),
    (x->-x[1], [randn(1) for i=1:3]),
    (x->x[1]^2, [randn(1) for i=1:3]),
    (x->x[1]^2.5, [rand(1) for i=1:3]),
    (x->1/x[1], [randn(1) for i=1:3]),
    (x->1/x[1], [randn(1) for i=1:3]),
    (x->sin(x[1]), [randn(1) for i=1:3]),
    (x->cos(x[1]), [randn(1) for i=1:3]),
    (x->exp(x[1]), [randn(1) for i=1:3]),
    (x->x[1]+x[2], [randn(2) for i=1:3]),
    (x->x[1]-x[2], [randn(2) for i=1:3]),
    (x->x[1]*x[2], [randn(2) for i=1:3]),
    (x->x[1]/x[2], [randn(2) for i=1:3]),
    (x->sin(x[1]/x[2]), [randn(2) for i=1:3]),
    (x->cos(x[1]-x[2])^2, [randn(2) for i=1:3]),
    (x->exp(2^x[1]+x[2]^3), [randn(2) for i=1:3]),
    (x->x[1]/log(x[2]^2+9.), [randn(2) for i=1:3]),
]
    @test compare(f,func(f(Source())), points)
end
