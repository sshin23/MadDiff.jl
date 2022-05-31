using Distributed, Plots, BenchmarkProfiles, DelimitedFiles; pgfplotsx()

const NWORKERS = 10

if nworkers() != NWORKERS
    addprocs(NWORKERS - nworkers(),exeflags="--project=.")
end

@everywhere begin

    using MadDiff, JuMP, Ipopt, AmplNLWriter, Ipopt_jll

    excludes = [
        "brainpc0.nl","brainpc1.nl","brainpc2.nl","brainpc3.nl","brainpc4.nl","brainpc5.nl","brainpc6.nl","brainpc7.nl","brainpc8.nl",
        "cresc132.nl",
        "drcav1lq.nl",
        "drcav2lq.nl","drcav3lq.nl","drcavty1.nl","drcavty2.nl","drcavty3.nl","sensors.nl"
    ]
    
    cases = readdir(ENV["CUTE_DIR"])
    setdiff!(cases, excludes)

    function replace_ad_backend(m,ad)
        m.moi_backend.model_cache.modattr[MOI.NLPBlock()] =
            MOI.NLPBlockData(MOI.Nonlinear.Evaluator(
                m.moi_backend.model_cache.modattr[MOI.NLPBlock()].evaluator.model,
                ad,
                Vector{MOI.VariableIndex}()));
    end

    function get_t3(model)
        return model.moi_backend.optimizer.model.nlp_data.evaluator.eval_objective_timer +
            model.moi_backend.optimizer.model.nlp_data.evaluator.eval_constraint_timer +
            model.moi_backend.optimizer.model.nlp_data.evaluator.eval_objective_gradient_timer +
            model.moi_backend.optimizer.model.nlp_data.evaluator.eval_constraint_jacobian_timer +
            model.moi_backend.optimizer.model.nlp_data.evaluator.eval_hessian_lagrangian_timer
    end


    function benchmark_cute(case)
        println("$case")
        
        (t1_jump,t2_jump,t3_jump) =  try 
            m = read_from_file(joinpath(ENV["CUTE_DIR"],case))
            set_optimizer(m,Ipopt.Optimizer)
            set_optimizer_attribute(m,"linear_solver","ma27")
            optimize!(m)

            m = read_from_file(joinpath(ENV["CUTE_DIR"],case))
            set_optimizer(m,Ipopt.Optimizer)
            set_optimizer_attribute(m,"linear_solver","ma27")
            GC.enable(false)
            t1 = @elapsed begin
                optimize!(m)
            end
            GC.enable(true)
            termination_status(m) == MOI.LOCALLY_SOLVED  || error()
            
            (t1, solve_time(m), get_t3(m))
        catch e
            (Inf, Inf, Inf)
        end
        
        (t1_maddiff,t2_maddiff,t3_maddiff) = try
            m = read_from_file(joinpath(ENV["CUTE_DIR"],case))
            set_optimizer(m,Ipopt.Optimizer)
            set_optimizer_attribute(m,"linear_solver","ma27")
            replace_ad_backend(m, MadDiffMOI.MadDiffAD())        
            optimize!(m);

            m = read_from_file(joinpath(ENV["CUTE_DIR"],case))
            set_optimizer(m,Ipopt.Optimizer)
            set_optimizer_attribute(m,"linear_solver","ma27")
            replace_ad_backend(m, MadDiffMOI.MadDiffAD())
            GC.enable(false)
            t1 = @elapsed begin
                optimize!(m)
            end
            GC.enable(true)
            termination_status(m) == MOI.LOCALLY_SOLVED || error()

            (t1, solve_time(m), get_t3(m))
            
        catch e
            (Inf, Inf, Inf)
        end


        (t1_ampl,t2_ampl,t3_ampl) = try
            m = read_from_file(joinpath(ENV["CUTE_DIR"],case))
            set_optimizer(m,() -> AmplNLWriter.Optimizer(Ipopt_jll.amplexe));
            set_optimizer_attribute(m,"linear_solver","ma27")
            optimize!(m) # force compile

            m = read_from_file(joinpath(ENV["CUTE_DIR"],case))
            set_optimizer(m,() -> AmplNLWriter.Optimizer(Ipopt_jll.amplexe));
            set_optimizer_attribute(m,"linear_solver","ma27")
            GC.enable(false)
            t1 = @elapsed begin
                optimize!(m)
            end
            GC.enable(true)
            
            termination_status(m) == MOI.LOCALLY_SOLVED  || error()

            (t1, solve_time(m), NaN)
        catch e
            (Inf, Inf, Inf)
        end
        return (
            t1_jump, t2_jump, t3_jump,
            t1_maddiff, t2_maddiff, t3_maddiff,
            t1_ampl, t2_ampl
        )
    end
end

results = pmap(benchmark_cute, cases[1:10])

t1s_jump = Float64[]
t2s_jump = Float64[]
t3s_jump = Float64[]
t2s_maddiff = Float64[]
t1s_maddiff = Float64[]
t3s_maddiff = Float64[]
t1s_ampl = Float64[]
t2s_ampl = Float64[]
# t3s_ampl = Float64[]

for (t1_jump, t2_jump, t3_jump, t1_maddiff, t2_maddiff, t3_maddiff, t1_ampl, t2_ampl) in results
    push!(t1s_jump, t1_jump)
    push!(t1s_maddiff, t1_maddiff)
    push!(t1s_ampl, t1_ampl)
    push!(t2s_jump, t2_jump)
    push!(t2s_maddiff, t2_maddiff)
    push!(t2s_ampl, t2_ampl)
    push!(t3s_jump, t3_jump)
    push!(t3s_maddiff, t3_maddiff)
end


writedlm("t1s_jump.csv", t1s_jump, ',')
writedlm("t2s_jump.csv", t2s_jump, ',')
writedlm("t3s_jump.csv", t3s_jump, ',')

writedlm("t1s_maddiff.csv", t1s_maddiff, ',')
writedlm("t2s_maddiff.csv", t2s_maddiff, ',')
writedlm("t3s_maddiff.csv", t3s_maddiff, ',')

writedlm("t1s_ampl.csv", t1s_ampl, ',')
writedlm("t2s_ampl.csv", t2s_ampl, ',')
# writedlm("t3s_ampl.csv", t3s_ampl, ',')

p1 = performance_profile(PlotsBackend(), Matrix{Float64}([t1s_jump t1s_maddiff t1s_ampl]), ["JuMP", "MadDiff", "Ampl"], title="optimize! time", logscale =false, framestyle=:box, legend=:bottomright)
p2 = performance_profile(PlotsBackend(), Matrix{Float64}([t2s_jump t2s_maddiff t2s_ampl]), ["JuMP", "MadDiff", "Ampl"], title="solver time", logscale =false, framestyle=:box, legend=:bottomright)
p3 = performance_profile(PlotsBackend(), Matrix{Float64}([t3s_jump t3s_maddiff]), ["JuMP", "MadDiff"], title="AD time", logscale =false, framestyle=:box, legend=:bottomright)

savefig(p1,"opt-time.pdf")
savefig(p2,"solver-time.pdf")
savefig(p3,"ad-time.pdf")
