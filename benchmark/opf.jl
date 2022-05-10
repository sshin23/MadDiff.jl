using MadDiff, PowerModels, JuMP, Ipopt, Plots
pgfplotsx()
PowerModels.silence()

# cases = [
#     "pglib_opf_case24_ieee_rts.m", "pglib_opf_case4619_goc.m", "pglib_opf_case2736sp_k.m"
# ]

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



function get_t3(model)
    return model.moi_backend.optimizer.model.nlp_data.evaluator.eval_objective_timer +
        model.moi_backend.optimizer.model.nlp_data.evaluator.eval_constraint_timer +
        model.moi_backend.optimizer.model.nlp_data.evaluator.eval_objective_gradient_timer +
        model.moi_backend.optimizer.model.nlp_data.evaluator.eval_constraint_jacobian_timer +
        model.moi_backend.optimizer.model.nlp_data.evaluator.eval_hessian_lagrangian_timer
end

t1s_jump = []
t1s_maddiff = []
t2s_jump = []
t2s_maddiff = []
t3s_jump = []
t3s_maddiff = []

for case in cases
    m = instantiate_model(
        joinpath("$(ENV["PGLIB_PATH"])",case),
        ACPPowerModel,
        PowerModels.build_opf
    ).model

    set_optimizer(m,Ipopt.Optimizer)
    set_optimizer_attribute(m,"linear_solver","ma27")
    
    t1_jump = @elapsed begin
        optimize!(m)
    end
    t2_jump = solve_time(m)
    t3_jump = get_t3(m)

    t1_maddiff = @elapsed begin
        optimize!(m,differentiation_backend = MadDiffAutomaticDifferentiation())
    end
    t2_maddiff = solve_time(m)
    t3_maddiff = get_t3(m)

    push!(t1s_jump, t1_jump)
    push!(t1s_maddiff, t1_maddiff)
    push!(t2s_jump, t2_jump)
    push!(t2s_maddiff, t2_maddiff)
    push!(t3s_jump, t3_jump)
    push!(t3s_maddiff, t3_maddiff)
end


p1 = performance_profile(plt, [t1s_jump t1s_maddiff], ["JuMP", "MadDiff"], title="optimize! time", logscale =false)
p2 = performance_profile(plt, [t2s_jump t2s_maddiff], ["JuMP", "MadDiff"], title="solver time", logscale =false)
p3 = performance_profile(plt, [t3s_jump t3s_maddiff], ["JuMP", "MadDiff"], title="AD time", logscale =false)

savefiv(p1,"opt-time.pdf")
savefiv(p2,"solver-time.pdf")
savefiv(p3,"ad-time.pdf")
