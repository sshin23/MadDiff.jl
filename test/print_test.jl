# We only do a rought test of printing

@test begin

    m = SimpleNLModels.Model()
    u = [variable(m;name = "u[$i]") for i=1:10]
    v = [parameter(m;name = "v[$i]") for i=1:10]
    
    println(u[2])
    println(v[3])
    
    for (x,p) in [(Variable("xx"),Parameter("pp")), (u,v)]

        println(sum(1 +beta(sin(x[1]+erf(p[2])/p[1])^2,cos(p[1]*x[4])/x[5]^2-x[3])+sin(x[i]) for i=1:7))
        println(1 - p[1] - sum(-sin(sin(sin(cos(beta(1,x[i])+p[2])-1.)+4)^9)/abs2(x[3]) for i=1:7))
        
    end
    
    true
end
