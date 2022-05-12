nlp_test = Dict()

nlp_test[1] = function (optimize!;opt...)

    m = MadDiffModel(;opt...)

    x = variable(m;lb=2,ub=4)
    objective(m,(x-1)^2)

    instantiate!(m)
    optimize!(m)
    
    return compare(value(x), 2.)
end

nlp_test[2] = function (optimize!;opt...)

    m = MadDiffModel(;opt...)

    x = variable(m)
    constraint(m,x+1<= -1)
    constraint(m,x+1>= -1)
    constraint(m,-2<= x+1)
    constraint(m,0>=x+1)

    instantiate!(m)
    optimize!(m)
    
    return compare(value(x), -2)
end

nlp_test[3] = function (optimize!;opt...)

    m = MadDiffModel(;opt...)

    x = [variable(m,name="s[$i]",start = .1) for i=1:3]
    
    c = constraint(m,0 ==x[1]+sin(x[2])-x[3]/2 + 1. )
    objective(m,x[2]^2 + x[2]^2 + 1.)

    instantiate!(m)
    optimize!(m)
    
    return compare(value.(x), [-0.74,0,.52]) && compare(dual(c),0.)
end

nlp_test[4] = function (optimize!;opt...)
    
    m = MadDiffModel(;opt...)

    x = [variable(m;start=mod(i,2)==1 ? -1.2 : 1.) for i=1:10];
    p = [parameter(m,2) for i=1:10];
    setvalue.(p,2)
    
    objective(m,sum(100(x[i-1]^2-x[i])^2+(x[i-1]-1)^2 for i=2:10))

    c = [ constraint(m,3x[i+1]^3+p[i]*x[i+2]-5+sin(x[i+1]-x[i+2])sin(x[i+1]+x[i+2])+4x[i+1]-x[i]exp(x[i]-x[i+1])-3 == 0) for i=1:8 ];

    instantiate!(m)
    optimize!(m)
    return compare(value.(x), [-0.95055636, 0.91390082, 0.98909052, 0.99855924, 0.99980874, 0.99997459, 0.99999662, 0.99999955, 0.99999994, 0.99999993])
end



m = MadDiffModel()

x = variable(m)
y = variable(m)
setvalue(x,3)
set_lower_bound(x,-1)
set_upper_bound(x,1)
@test MadDiffModels.index(x) == 1
@test value(x) == 3 
@test lower_bound(x) == -1 
@test upper_bound(x) == 1

c = constraint(m,x==1.)
constraint(m,x>=1.)
constraint(m,x<=1.)

constraint(m,1. ==y)
constraint(m,1. <=y)
constraint(m,-1. >=y)

constraint(m,x==y)
constraint(m,x-4. <=y)
constraint(m,x+4. >=y)

set_lower_bound(c,-1)
set_upper_bound(c,1)
@test lower_bound(c) == -1 
@test upper_bound(c) == 1

p = parameter(m,1.)
setvalue(p,4)
@test value(p) == 4

objective(m,sum([p,x,y]))

m[:test] = "test"

@test objective_value(m) == 7.
@test num_variables(m) == 2
@test num_constraints(m) == 9
@test m[:test] == "test"

for (optimizer,opt) in [(ipopt,[:print_level=>0]),
                        (madnlp,[:print_level=>MadNLP.ERROR])]
    for f in values(nlp_test)
        @test begin
            f(optimizer;opt...)
            true
        end
    end
end

# Test printing
@test begin
    println("Test printing")
    x = Variable()
    p = Parameter()
    m = MadDiffModel()
    u = [variable(m) for i=1:10]
    v = [parameter(m) for i=1:10]
    c = constraint(m,u[1]+u[2]==0)
    
    show(stdout, MIME"text/plain"(),x)
    show(stdout, MIME"text/plain"(),p)
    show(stdout, MIME"text/plain"(),x[1])
    show(stdout, MIME"text/plain"(),p[1])
    
    show(stdout, MIME"text/plain"(),2 + 7*sum(isodd(i) ? p[i]-(x[2]+((3+x[2])*((x[1]+p[2])+(Constant(2.)+sin(x[i]))*2) + (3+x[i])/p[1] - (x[1]+x[i])) + (x[2]+2)^(p[2]+1)) :
                                              1 +beta(sin(x[1]^x[2])+(erf(p[2])+sin(x[2])/p[1])^2,(cos(p[1]*x[4])+x[1])/x[5]^2-x[3]+sin(x[i])) + beta(p[2]*x[1],3) for i=1:7))
    show(stdout, MIME"text/plain"(),1 - (1/p[1])^2 - sum(-sin(sin(sin(cos(beta(1,x[i]/x[i+1])+p[2])-1.)+4)^9)/abs2(x[3]) for i=1:7))
    show(stdout, MIME"text/plain"(),sum(1 - (p[i]+x[i])^(p[2]+x[1]) + (p[i]+x[i])/(p[2]+x[1]) for i=1:7))
    show(stdout, MIME"text/plain"(),1 - p[2]^(p[1]+x[1]))
    show(stdout, MIME"text/plain"(),1 - p[2]/(p[1]+x[1]))
    show(stdout, MIME"text/plain"(),1 - p[2]*(p[1]+x[1]))
    show(stdout, MIME"text/plain"(),1 - (p[2]+x[3])+(p[1]+x[1]))
    
    true
end
