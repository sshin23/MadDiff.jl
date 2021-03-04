using SimpleNLModels, Ipopt, MadNLP, JuMP, AmplNLWriter, StatsPlots

include("benchmark_include.jl")

result = []

for (N,nlm,jump,name) in [
    (10000,nlm_luksan_vlcek_501,jump_luksan_vlcek_501,"luksanvlcek"),
    # (10000,nlm_ocp,jump_ocp,"hehnandrea")
]
    for i=1:2
        t1 = @elapsed begin
            @time begin
                m = nlm(;N=N,optimizer=Ipopt.Optimizer,output_file="output/$name-simplenlmodels.out")
                instantiate!(m)
            end
            @time optimize!(m)
        end
        # t2,t3 = parse_ipopt_output("output/$name-simplenlmodels.out")

        # m = jump(;N=N,optimizer=() -> AmplNLWriter.Optimizer("ipopt",["output_file=output/$name-ampl.out"]))
        # optimize!(m)
        # t4 = m.moi_backend.optimizer.model.ext[:MPBModel].inner.solve_time
        # t5,t6 = parse_ipopt_output("output/$name-ampl.out")

        # t7 = @elapsed begin
        #     m = jump(;N=N,optimizer=()->Ipopt.Optimizer(output_file="output/$name-jump.out"))
        #     optimize!(m)
        #     t8,t9 = parse_ipopt_output("output/$name-jump.out")
        # end

        # i==2 && push!(result,[[t1,t4,t7],[t2,t5,t8],[t3,t6,t9]])
    end
end

# ticklabel = ["SimpleNLModels","AMPL","JuMP"]
# barlabel = ["Model Creation + Solution Time (Wall)" "Solution Time (CPU)" "NLP Function Evaluation Time (CPU)"]

# plt1 = groupedbar([result[1][1] result[1][2] result[1][3]], title="Luksan & Vlcek (1999)", bar_position = :dodge,
#                   bar_width=0.7, xticks=(1:3, ticklabel), label=barlabel, ylabel="Time (sec)", framestyle=:box, ylim = (0,1.1*maximum(vcat(result[1]...))),dpi = 300)
# plt2 = groupedbar([result[2][1] result[2][2] result[2][3]], title="Hehn & Andrea (2011, modified)", bar_position = :dodge,
#                   bar_width=0.7, xticks=(1:3, ticklabel), label=barlabel, ylabel="Time (sec)", framestyle=:box, ylim = (0,1.1*maximum(vcat(result[2]...))),dpi = 300)
# savefig(plt1,"output/luksanvlcek.png")
# savefig(plt2,"output/hehnandrea.png")



