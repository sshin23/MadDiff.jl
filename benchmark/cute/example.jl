using MadDiff, JuMP, Ipopt, DelimitedFiles

const CUTE_DIR = "/home/sshin/git/nl-benchmark/cute/text/"

cases = readdir(CUTE_DIR)

# dual2.nl: summation compilation issue 

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

t2s_jump = Float64[]
t3s_jump = Float64[]
t2s_maddiff = Float64[]
t3s_maddiff = Float64[]

cnt = 0
for case in cases
    global cnt
    cnt += 1
    println("$cnt: $case")

    if cnt < current
        continue
    end
    
    m = try 
        m = read_from_file(joinpath(CUTE_DIR,case))
        set_optimizer(m,Ipopt.Optimizer)
        set_optimizer_attribute(m,"linear_solver","ma27")
        optimize!(m)

        m = read_from_file(joinpath(CUTE_DIR,case))
        set_optimizer(m,Ipopt.Optimizer)
        set_optimizer_attribute(m,"linear_solver","ma27")
        GC.enable(false)
        optimize!(m)
        GC.enable(true)
        solved = termination_status(m) == MOI.LOCALLY_SOLVED
        push!(t2s_jump, solved ? solve_time(m) : Inf)
        push!(t3s_jump, solved ? get_t3(m) : Inf)
        
        m
    catch e
        continue
    end
    
    try
        m = read_from_file(joinpath(CUTE_DIR,case))
        set_optimizer(m,Ipopt.Optimizer)
        set_optimizer_attribute(m,"linear_solver","ma27")
        replace_ad_backend(m, MadDiffMOI.MadDiffAD())        
        optimize!(m);

        m = read_from_file(joinpath(CUTE_DIR,case))
        set_optimizer(m,Ipopt.Optimizer)
        set_optimizer_attribute(m,"linear_solver","ma27")
        replace_ad_backend(m, MadDiffMOI.MadDiffAD())
        GC.enable(false)
        optimize!(m)
        GC.enable(true)
        solved = termination_status(m) == MOI.LOCALLY_SOLVED
        push!(t2s_maddiff, solved ? solve_time(m) : Inf)
        push!(t3s_maddiff, solved ? get_t3(m) : Inf)
        
    catch e
        push!(t2s_maddiff, Inf)
        push!(t3s_maddiff, Inf)
        
        throw(e)
    end
end



# using Distributed, Plots, BenchmarkProfiles

# const NWORKERS = 10

# if nworkers() != NWORKERS
#     addprocs(NWORKERS - nworkers(),exeflags="--project=.")
# end

# @everywhere begin
#     using MadDiff, PowerModels, JuMP, Ipopt, AmplNLWriter, Ipopt_jll, DelimitedFiles
#     PowerModels.silence()

#     cases = [
#         "pglib_opf_case24_ieee_rts.m", "pglib_opf_case4619_goc.m", "pglib_opf_case2736sp_k.m",
#         "pglib_opf_case4661_sdet.m", "pglib_opf_case2737sop_k.m", "pglib_opf_case4837_goc.m",
#         "pglib_opf_case2742_goc.m", "pglib_opf_case4917_goc.m", "pglib_opf_case2746wop_k.m",
#         "pglib_opf_case500_goc.m",  "pglib_opf_case2746wp_k.m", "pglib_opf_case57_ieee.m",
#         "pglib_opf_case10000_goc.m", "pglib_opf_case2848_rte.m", "pglib_opf_case588_sdet.m",
#         "pglib_opf_case10480_goc.m", "pglib_opf_case2853_sdet.m", "pglib_opf_case5_pjm.m",
#         "pglib_opf_case118_ieee.m", "pglib_opf_case2868_rte.m", "pglib_opf_case60_c.m",
#         "pglib_opf_case1354_pegase.m", "pglib_opf_case2869_pegase.m", "pglib_opf_case6468_rte.m",
#         "pglib_opf_case13659_pegase.m", "pglib_opf_case30000_goc.m", "pglib_opf_case6470_rte.m",
#         "pglib_opf_case14_ieee.m", "pglib_opf_case300_ieee.m", "pglib_opf_case6495_rte.m",
#         "pglib_opf_case162_ieee_dtc.m", "pglib_opf_case3012wp_k.m", "pglib_opf_case6515_rte.m",
#         "pglib_opf_case179_goc.m", "pglib_opf_case3022_goc.m", "pglib_opf_case73_ieee_rts.m",
#         "pglib_opf_case1888_rte.m", "pglib_opf_case30_as.m", "pglib_opf_case793_goc.m",
#         "pglib_opf_case19402_goc.m", "pglib_opf_case30_ieee.m", "pglib_opf_case8387_pegase.m",
#         "pglib_opf_case1951_rte.m", "pglib_opf_case3120sp_k.m", "pglib_opf_case89_pegase.m",
#         "pglib_opf_case2000_goc.m", "pglib_opf_case3375wp_k.m", "pglib_opf_case9241_pegase.m",
#         "pglib_opf_case200_activ.m", "pglib_opf_case3970_goc.m", "pglib_opf_case9591_goc.m",
#         "pglib_opf_case2312_goc.m", "pglib_opf_case39_epri.m", "pglib_opf_case2383wp_k.m",
#         "pglib_opf_case3_lmbd.m", "pglib_opf_case240_pserc.m", "pglib_opf_case4020_goc.m",
#         "pglib_opf_case24464_goc.m", "pglib_opf_case4601_goc.m"
#     ]





#     function benchmark_opf(case)
#         m = instantiate_model(
#             joinpath("$(ENV["PGLIB_PATH"])",case),
#             ACPPowerModel,
#             PowerModels.build_opf
#         ).model

#         # JUMP
#         set_optimizer(m,Ipopt.Optimizer)
#         set_optimizer_attribute(m,"linear_solver","ma27")
        
#         optimize!(m) # force compile
#         GC.enable(false)
#         t1_jump = @elapsed begin
#             optimize!(m)
#         end
#         GC.enable(true)
#         t2_jump = solve_time(m)
#         t3_jump = get_t3(m)

#         # MadDiff
#         set_optimizer(m,Ipopt.Optimizer)
#         set_optimizer_attribute(m,"linear_solver","ma27")
#         optimize!(m,differentiation_backend = MadDiffAD()) # force compile
#         GC.enable(false)
#         t1_maddiff = @elapsed begin
#             optimize!(m,differentiation_backend = MadDiffAD())
#         end
#         GC.enable(true)
#         t2_maddiff = solve_time(m)
#         t3_maddiff = get_t3(m)

#         # Ampl
#         set_optimizer(m,() -> AmplNLWriter.Optimizer(Ipopt_jll.amplexe));
#         set_optimizer_attribute(m,"linear_solver","ma27")

#         optimize!(m) # force compile
#         GC.enable(false)
#         t1_ampl = @elapsed begin
#             optimize!(m,differentiation_backend = MOI.Nonlinear.ExprGraphOnly())
#         end
#         GC.enable(true)
#         t2_ampl = solve_time(m)
#         # t3_ampl = get_t3(m)

        
#         # push!(t3s_ampl, t3_ampl)

#         return (
#             t1_jump, t2_jump, t3_jump,
#             t1_maddiff, t2_maddiff, t3_maddiff,
#             t1_ampl, t2_ampl
#         )
#     end

# end

# results = pmap(benchmark_opf, cases)

# t1s_jump = Float64[]
# t2s_jump = Float64[]
# t3s_jump = Float64[]
# t2s_maddiff = Float64[]
# t1s_maddiff = Float64[]
# t3s_maddiff = Float64[]
# t1s_ampl = Float64[]
# t2s_ampl = Float64[]
# # t3s_ampl = Float64[]

# for (t1_jump, t2_jump, t3_jump, t1_maddiff, t2_maddiff, t3_maddiff, t1_ampl, t2_ampl) in results
#     push!(t1s_jump, t1_jump)
#     push!(t1s_maddiff, t1_maddiff)
#     push!(t1s_ampl, t1_ampl)
#     push!(t2s_jump, t2_jump)
#     push!(t2s_maddiff, t2_maddiff)
#     push!(t2s_ampl, t2_ampl)
#     push!(t3s_jump, t3_jump)
#     push!(t3s_maddiff, t3_maddiff)
# end


# writedlm("t1s_jump.csv", t1s_jump, ',')
# writedlm("t2s_jump.csv", t2s_jump, ',')
# writedlm("t3s_jump.csv", t3s_jump, ',')

# writedlm("t1s_maddiff.csv", t1s_maddiff, ',')
# writedlm("t2s_maddiff.csv", t2s_maddiff, ',')
# writedlm("t3s_maddiff.csv", t3s_maddiff, ',')

# writedlm("t1s_ampl.csv", t1s_ampl, ',')
# writedlm("t2s_ampl.csv", t2s_ampl, ',')
# # writedlm("t3s_ampl.csv", t3s_ampl, ',')

# t1s_symb = readdlm("t1s_symb.csv",Float64)
# t2s_symb = readdlm("t2s_symb.csv",Float64)
# t3s_symb = readdlm("t3s_symb.csv",Float64)

# p1 = performance_profile(PlotsBackend(), Matrix{Float64}([t1s_jump t1s_maddiff t1s_symb t1s_ampl]), ["JuMP", "MadDiff", "SymbolicAD", "Ampl"], title="optimize! time", logscale =false, framestyle=:box, legend=:bottomright)
# p2 = performance_profile(PlotsBackend(), Matrix{Float64}([t2s_jump t2s_maddiff t2s_symb t2s_ampl]), ["JuMP", "MadDiff", "SymbolicAD", "Ampl"], title="solver time", logscale =false, framestyle=:box, legend=:bottomright)
# p3 = performance_profile(PlotsBackend(), Matrix{Float64}([t3s_jump t3s_maddiff t3s_symb]), ["JuMP", "MadDiff", "SymbolicAD"], title="AD time", logscale =false, framestyle=:box, legend=:bottomright)

# savefig(p1,"opt-time.pdf")
# savefig(p2,"solver-time.pdf")
# savefig(p3,"ad-time.pdf")
