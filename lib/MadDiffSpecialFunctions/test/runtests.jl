using MadDiffSpecialFunctions, SpecialFunctions, Test, MadDiffTests
import Random: seed!

seed!(0)

@testset "MadDiff test" begin
    xs = [rand(10) for i=1:3]
    ps = [rand(10) for i=1:3]
    
    y1 = zeros(20);
    y2 = zeros(20);
    z1 = zeros(10,10);
    z2 = zeros(10,10);



    test_list_1 = [
        ("function-test-1-1",(x,p)->beta(erf(x[1]/x[2]/p[3])+p[3]*x[2],erf(x[9])^2)), 
        ("function-test-1-3",(x,p)->beta(cos(log(inv(inv(x[1])))),erfc(tanh(0*x[1])))), 
        ("function-test-1-6",(x,p)->beta(2*logbeta(x[1],x[5]),beta(x[2],x[3]))), 
        ("function-test-1-7",(x,p)->besselj0(exp(erf(-x[1])))), 
        ("function-test-1-8",(x,p)->erfc((x[1]^2/x[2])^x[9]/x[10])), 
        ("function-test-1-9",(x,p)->erfc(x[1])^erf(2.5x[2])), 
        ("function-test-1-12",(x,p)->sum(isodd(i) ? 4+x[i] : sin(x[1]) for i=1:4)*erf(x[1]/p[1])/p[10]/p[9]*p[5]), 
        ("function-test-1-14",(x,p)->airyai(exp(x[1]+x[2]*p[2]^8))), 
        ("function-test-1-21",(x,p)->beta(beta(tan(beta(x[1],1)+p[2]),cos(sin(x[2]))),x[3])), 
        ("function-test-1-23",(x,p)->sum(1 +beta(sin(x[1]+erf(p[2])/p[1])^2,cos(p[1]*x[4])/x[5]^2-x[3])+sin(x[i]) for i=1:7)),
        ("function-test-1-24",(x,p)->beta(cos(beta(beta(x[1]^9,x[2]),x[2]*x[3])),sin(x[2]*x[3]/p[2])/p[1]))
    ]

    test_list_2 = [
        ("function-test-2-2",(y,x,p)->begin
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

    for f in MadDiffSpecialFunctions._UNIVARIATE_FUNCTIONS
        forig = (x,p)->f(x[1])
        name = "basic-function-test-$f"
        MadDiffTests.test_function(forig,name, xs, ps,y1,y2,z1,z2)
    end

    for (name,forig) in test_list_1
        MadDiffTests.test_function(forig,name,xs,ps,y1,y2,z1,z2)
    end

    for (name,forig) in test_list_2
        MadDiffTests.test_field(forig,name,xs,ps,y1,y2,z1,z2)
    end
end 
