using Test, SimpleNLModels, MadNLP
import Random: seed!
const eps = 1e-6
seed!(0)

compare(p1,p2) = maximum(abs.(p1 .- p2)) <= eps

function compare(f1,f2,points)
    for point in points
        compare(f1(point),f2(point)) || return false
    end
    return true
end

@testset "SimpleNLModels test" begin
    @testset "Function Test" begin
        include("function_test.jl")
    end
    @testset "NLP Test" begin
        include("nlp_test.jl")
    end
    @testset "Ipopt Test" begin
        include("ipopt_test.jl")
    end
    @testset "Madnlp Test" begin
        include("madnlp_test.jl")
    end
end 
