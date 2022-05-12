using Distributed

const NWORKERS = 10

if nworkers() != NWORKERS
    addprocs(NWORKERS - nworkers(),exeflags="--project=.")
end

cases = [
    "pglib_opf_case24_ieee_rts.m", "pglib_opf_case4619_goc.m", "pglib_opf_case2736sp_k.m",
    "pglib_opf_case4661_sdet.m", "pglib_opf_case2737sop_k.m", "pglib_opf_case4837_goc.m",
    "pglib_opf_case2742_goc.m", "pglib_opf_case4917_goc.m", "pglib_opf_case2746wop_k.m",
    "pglib_opf_case500_goc.m",  "pglib_opf_case2746wp_k.m", "pglib_opf_case57_ieee.m",
    "pglib_opf_case10000_goc.m", "pglib_opf_case2848_rte.m", "pglib_opf_case588_sdet.m",
    "pglib_opf_case10480_goc.m", "pglib_opf_case2853_sdet.m", "pglib_opf_case5_pjm.m",
    "pglib_opf_case118_ieee.m", "pglib_opf_case2868_rte.m", "pglib_opf_case60_c.m",
    "pglib_opf_case1354_pegase.m", "pglib_opf_case2869_pegase.m", "pglib_opf_case6468_rte.m",
    "pglib_opf_case13659_pegase.m", "pglib_opf_case30000_goc.m", "pglib_opf_case6470_rte.m",
    "pglib_opf_case14_ieee.m", "pglib_opf_case300_ieee.m", "pglib_opf_case6495_rte.m",
    "pglib_opf_case162_ieee_dtc.m", "pglib_opf_case3012wp_k.m", "pglib_opf_case6515_rte.m",
    "pglib_opf_case179_goc.m", "pglib_opf_case3022_goc.m", "pglib_opf_case73_ieee_rts.m",
    "pglib_opf_case1888_rte.m", "pglib_opf_case30_as.m", "pglib_opf_case793_goc.m",
    "pglib_opf_case19402_goc.m", "pglib_opf_case30_ieee.m", "pglib_opf_case8387_pegase.m",
    "pglib_opf_case1951_rte.m", "pglib_opf_case3120sp_k.m", "pglib_opf_case89_pegase.m",
    "pglib_opf_case2000_goc.m", "pglib_opf_case3375wp_k.m", "pglib_opf_case9241_pegase.m",
    "pglib_opf_case200_activ.m", "pglib_opf_case3970_goc.m", "pglib_opf_case9591_goc.m",
    "pglib_opf_case2312_goc.m", "pglib_opf_case39_epri.m", "pglib_opf_case2383wp_k.m",
    "pglib_opf_case3_lmbd.m", "pglib_opf_case240_pserc.m", "pglib_opf_case4020_goc.m",
    "pglib_opf_case24464_goc.m", "pglib_opf_case4601_goc.m"
]

@everywhere begin
    using SymbolicAD, PowerModels, JuMP, Ipopt, DelimitedFiles
    PowerModels.silence()

    function get_t3(model)
        return model.moi_backend.optimizer.model.nlp_data.evaluator.eval_objective_timer +
            model.moi_backend.optimizer.model.nlp_data.evaluator.eval_constraint_timer +
            model.moi_backend.optimizer.model.nlp_data.evaluator.eval_objective_gradient_timer +
            model.moi_backend.optimizer.model.nlp_data.evaluator.eval_constraint_jacobian_timer +
            model.moi_backend.optimizer.model.nlp_data.evaluator.eval_hessian_lagrangian_timer
    end
    
    function benchmark_opf(case)
        m = instantiate_model(
            joinpath("$(ENV["PGLIB_PATH"])",case),
            ACPPowerModel,
            PowerModels.build_opf
        ).model

        set_optimizer(m,Ipopt.Optimizer)
        set_optimizer_attribute(m,"linear_solver","ma27")
        set_optimize_hook(m, SymbolicAD.optimize_hook)
        
        optimize!(m) # force compile
        GC.enable(false)
        t1_symb = @elapsed begin
            optimize!(m)
        end
        GC.enable(true)
        t2_symb = solve_time(m)
        t3_symb = get_t3(m)

        return t1_symb, t2_symb, t3_symb
    end
end
    

results = pmap(benchmark_opf,cases)

t1s_symb = []
t2s_symb = []
t3s_symb = []
for (t1_symb,t2_symb,t3_symb) in results
    push!(t1s_symb, t1_symb)
    push!(t2s_symb, t2_symb)
    push!(t3s_symb, t3_symb)
end

writedlm("../t1s_symb.csv", t1s_symb, ',')
writedlm("../t2s_symb.csv", t2s_symb, ',')
writedlm("../t3s_symb.csv", t3s_symb, ',')
