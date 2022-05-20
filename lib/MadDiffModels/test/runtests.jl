using Test, SpecialFunctions, NLPModelsIpopt, MadNLP
using MadDiffModels: MadDiffModel, variable, parameter, constraint, objective, value, setvalue, lower_bound, upper_bound, set_lower_bound, set_upper_bound, instantiate!, get_nvar, get_ncon
const eps = 1e-6

compare(p1,p2) = maximum(abs.(p1 .- p2)) <= eps

function compare(f1,f2,xpoints,ppoints)
    for i=1:length(xpoints)
        compare(f1(xpoints[i],ppoints[i]),f2(xpoints[i],ppoints[i])) || return false
    end
    return true
end

function compare(f1,f2,y1,y2,xpoints,ppoints)
    for i=1:length(xpoints)
        f1(y1,xpoints[i],ppoints[i])
        f2(y2,xpoints[i],ppoints[i])
        compare(y1,y2) || return false
    end
    return true
end

@testset "MadDiffModels test" begin
    @testset "Function Test" begin
        include("model_test.jl")
    end
end 
