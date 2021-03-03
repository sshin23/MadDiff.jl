nlp_test = Dict()

nlp_test[1] = function (optimizer;opt...)

    m = SimpleNLModels.Model(optimizer;opt...)

    x = variable(m;lb=2,ub=4)
    objective(m,(x-1)^2)

    optimize!(m)
    
    return compare(value(x), 2.)
end

nlp_test[2] = function (optimizer;opt...)

    m = SimpleNLModels.Model(optimizer;opt...)

    x = variable(m)
    constraint(m,x+1,lb=-1,ub=-1)

    optimize!(m)
    
    return compare(value(x), -2)
end

nlp_test[3] = function (optimizer;opt...)

    m = SimpleNLModels.Model(optimizer;opt...)

    x = [variable(m,name="s[$i]",start = .1) for i=1:3]
    
    constraint(m,x[1]+sin(x[2])-x[3]/2 + 1.)
    objective(m,x[2]^2)

    optimize!(m)
    
    return compare(value.(x), [-0.74,0,.52])
end

for (optimizer,opt) in [(IpoptOptimizer,[:print_level=>0]),
                        (MadNLPOptimizer,[:print_level=>MadNLP.ERROR])]
    for f in values(nlp_test)
        @test f(optimizer;opt...)
    end
end
