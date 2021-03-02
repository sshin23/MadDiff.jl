using Test, NonlinearModels
import Random: seed!
import LinearAlgebra: norm
seed!(0)

function compare(f1,f2,points)
    for point in points
        norm(f1(point) .- f2(point)) <= 1e-8 || return false
    end
    return true
end

@testset "MadNLP test" begin
    @testset "Function Test" begin
        include("function_test.jl")
    end
    @testset "Print Test" begin
        include("print_test.jl")
    end
    @testset "NLP Test" begin
        include("nlp_test.jl")
    end
end # @testset
