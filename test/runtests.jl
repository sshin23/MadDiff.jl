using Test, SimpleNL, MadNLP, Ipopt, ForwardDiff, SpecialFunctions
import Random: seed!
const eps = 1e-6
seed!(0)

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

@testset "SimpleNLModels test" begin
    @testset "basic" begin
        include("basic_test.jl")
    end
    @testset "Function Test" begin
        include("function_test.jl")
    end
    @testset "Print Test" begin
        include("print_test.jl")
    end
    @testset "NLP Test" begin
        include("nlp_test.jl")
    end
end 
