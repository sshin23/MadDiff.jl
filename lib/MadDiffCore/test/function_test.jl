y1 = zeros(20);
y2 = zeros(20);
z1 = zeros(10,10);
z2 = zeros(10,10);

xpoints1 = [rand(10) for i=1:3]
xpoints2 = [rand(10).+1 for i=1:3]
ppoints = [rand(10) for i=1:3]

function force_lower_tri(z)
    for i=1:size(z,1)
        for j=1:size(z,2)
            i<j && (z[i,j]=0)
        end
    end
end

get_gorig!(forig) = (y,x,p) -> (y.=0;ForwardDiff.gradient!(y,x->forig(x,p),x))
get_sgorig!(forig,gsparsity) = (y,x,p) ->  (y.=0;y[1:length(gsparsity)] .= ForwardDiff.gradient(x->forig(x,p),x)[gsparsity])
get_horig!(forig) = (z,x,p) -> (z.=0;ForwardDiff.hessian!(z,x->forig(x,p),x); force_lower_tri(z))
get_shorig!(forig,hsparsity) = (z,x,p) ->  (z.=0;hess=ForwardDiff.hessian(x->forig(x,p),x); z[1:length(hsparsity)] .= [hess[i,j] for (i,j) in hsparsity])
get_jorig!(forig) = (z,x,p) -> (z.=0;ForwardDiff.jacobian!(z,(y,x)->forig(y,x,p),zeros(size(z,1)),x))
get_sjorig!(forig,jsparsity) = (y,x,p) ->  (y.=0;n=length(y);jac=ForwardDiff.jacobian((y,x)->forig(y,x,p),zeros(n),x); y[1:length(jsparsity)] .= [jac[i,j] for (i,j) in jsparsity])

const test_list_1 = [
    ("function-test-1-1",(x,p)->beta(erf(x[1]/x[2]/p[3])+p[3]*x[2],erf(x[9])^2)), 
    ("function-test-1-2",(x,p)->0*x[1]), 
    ("function-test-1-3",(x,p)->beta(cos(log(inv(inv(x[1])))),erfc(tanh(0*x[1])))), 
    ("function-test-1-4",(x,p)->(0*x[1]^x[3]^p[1]+x[1])/x[9]/x[10]), 
    ("function-test-1-5",(x,p)->(x[1]+1.)^x[2]*log(x[3])/tanh(x[2])), 
    ("function-test-1-6",(x,p)->beta(2*logbeta(x[1],x[5]),beta(x[2],x[3]))), 
    ("function-test-1-7",(x,p)->besselj0(exp(erf(-x[1])))), 
    ("function-test-1-8",(x,p)->erfc((x[1]^2/x[2])^x[9]/x[10])), 
    ("function-test-1-9",(x,p)->erfc(x[1])^erf(2.5x[2])), 
    ("function-test-1-10",(x,p)->sin(1/x[1])), 
    ("function-test-1-11",(x,p)->sum(isodd(i) ? sin(x[i])/cos(x[i+1]) : csc(x[i+1]/p[2]+x[i-1]/p[3]) for i=2:9)), 
    ("function-test-1-12",(x,p)->sum(isodd(i) ? 4+x[i] : sin(x[1]) for i=1:4)*erf(x[1]/p[1])/p[10]/p[9]*p[5]), 
    ("function-test-1-13",(x,p)->exp(x[2])/cos(x[1])^2+sin(x[1]^2)), 
    ("function-test-1-14",(x,p)->airyai(exp(x[1]+x[2]*p[2]^8))), 
    ("function-test-1-15",(x,p)->prod(x[1]+x[2] for i=1:4)), 
    ("function-test-1-16",(x,p)->sin(x[9]inv(x[1])-x[8]inv(x[2]))), 
    ("function-test-1-17",(x,p)->-(x[2]-x[1])*(x[1]*x[2])*(x[1]*x[2])/sum(x[i]+p[j] for i=1:3 for j=4:8)), 
    ("function-test-1-18",(x,p)->sum(sum(isodd(i) ? x[1]/x[2] : cos(x[2]) for i =1:3) / sum(isodd(i) ? x[i]^x[j] : x[i]^2 for i=1:4) for j=1:3)), 
    ("function-test-1-19",(x,p)->x[1]/log(x[2]^2+9.)), 
    ("function-test-1-20",(x,p)->sum(cos(x[i]) for i=1:4)), 
    ("function-test-1-21",(x,p)->beta(beta(tan(beta(x[1],1)+p[2]),cos(sin(x[2]))),x[3])), 
    ("function-test-1-22",(x,p)->sum(f(x[1]) for f in [sin, cos, exp])),
    ("function-test-1-23",(x,p)->sum(1 +beta(sin(x[1]+erf(p[2])/p[1])^2,cos(p[1]*x[4])/x[5]^2-x[3])+sin(x[i]) for i=1:7)),
    ("function-test-1-24",(x,p)->beta(cos(beta(beta(x[1]^9,x[2]),x[2]*x[3])),sin(x[2]*x[3]/p[2])/p[1]))
]

const test_list_2 = [
    ("function-test-2-1",(y,x,p)->begin 
     y[1] = x[1]/log(x[2]^2+9.);
     y[2] = sin(1/x[1]);
     y[3] = sin(exp(x[1]+x[2]*p[2]^8));
     end),
    ("function-test-2-2",(y,x,p)->begin
     y[1] = (0*x[1]^x[3]^p[1]+x[1])/x[9]/x[10];
     y[2] = beta(2*logbeta(x[1],x[5]),beta(x[2],x[3]));
     y[3] = besselj0(exp(erf(-x[1]))); 
     y[4] = erfc((x[1]^2/x[2])^x[9]/x[10]);
     end),
    ("function-test-2-3",(y,x,p)-> begin
     y[1] = erfc(x[1])^erf(2.5x[2]);
     y[2] = prod(x[1]+x[2] for i=1:4); 
     y[3] = sin(x[9]inv(x[1])-x[8]inv(x[2]));
     y[4] = -(x[2]-x[1])*(x[1]*x[2])*(x[1]*x[2])/sum(x[i]+p[j] for i=1:3 for j=4:8);
     y[5] = sum(sum(x[1]/x[2] for i =1:3) / sum(x[i]^x[j] for i=1:4) for j=1:3)
     end),
 ]

function test_function_1(forig,name,xpoints,ppoints)
    ff = forig(MadDiffCore.Variable(),MadDiffCore.Parameter())
    f  = MadDiffCore.function_evaluator(ff)
    
    g! = MadDiffCore.gradient_evaluator(ff)
    sg!,gsparsity = MadDiffCore.sparse_gradient_evaluator(ff)
    
    h! = MadDiffCore.hessian_evaluator(ff)
    sh!,hsparsity = MadDiffCore.sparse_hessian_evaluator(ff)
    
    gorig!  = get_gorig!(forig)
    sgorig! = get_sgorig!(forig,gsparsity)
    horig!  = get_horig!(forig)
    shorig!  = get_shorig!(forig,hsparsity)

    @testset "$name" begin
        @test compare(forig,f,xpoints,ppoints)
        @test compare(gorig!,g!,y1,y2,xpoints,ppoints)
        @test compare(sgorig!,sg!,y1,y2,xpoints,ppoints)
        @test compare(horig!,h!,z1,z2,xpoints,ppoints)
        @test compare(shorig!,sh!,z1,z2,xpoints,ppoints)
    end
end

function test_function_2(forig,name,xpoints,ppoints)
    ff = MadDiffCore.Field()
    forig(ff,MadDiffCore.Variable(),MadDiffCore.Parameter())
    
    f  = MadDiffCore.field_evaluator(ff)
    
    j! = MadDiffCore.jacobian_evaluator(ff)
    sj!,jsparsity = MadDiffCore.sparse_jacobian_evaluator(ff)
    
    jorig!  = get_jorig!(forig)
    sjorig! = get_sjorig!(forig,jsparsity)

    @testset "$name" begin
        @test compare(forig,f,y1,y2,xpoints,ppoints)
        @test compare(jorig!,j!,z1,z2,xpoints,ppoints)
        @test compare(sjorig!,sj!,y1,y2,xpoints,ppoints)
    end
end

for (f,~,~) in MadDiffCore.f_nargs_1
    forig = @eval (x,p)->$f(x[1])
    name = "basic-function-test-$f"
    test_function_1(forig,name, f != :acoth ? xpoints1 : xpoints2, ppoints)
end

for (name,forig) in test_list_1
    test_function_1(forig,name,xpoints1, ppoints)
end

for (name,forig) in test_list_2
    test_function_2(forig,name,xpoints1, ppoints)
end


