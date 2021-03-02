using SimpleNLModels, Ipopt, JuMP, AmplNLWriter, StatsPlots

include("benchmark_include.jl")

result = []

for (N,nlm,jump) in [
    (10000,nlm_luksan_vlcek_501,jump_luksan_vlcek_501),
    (10000,nlm_ocp,jump_ocp)
]
    for i=1:2
        t1 = @elapsed begin
            m = nlm(;N=N,optimizer=SimpleNLModels.IpoptOptimizer,output_file="test.out")
            instantiate!(m)
            optimize!(m)
        end
        t2,t3 = parse_ipopt_output("test.out")

        m = jump(;N=N,optimizer=() -> AmplNLWriter.Optimizer("ipopt",["output_file=test.out"]))
        optimize!(m)
        t4 = m.moi_backend.optimizer.model.ext[:MPBModel].inner.solve_time
        t5,t6 = parse_ipopt_output("test.out")

        t7 = @elapsed begin
            m = jump(;N=N,optimizer=()->Ipopt.Optimizer(output_file="test.out"))
            optimize!(m)
            t8,t9 = parse_ipopt_output("test.out")
        end

        i==2 && push!(result,[[t1,t4,t7],[t2,t5,t8],[t3,t6,t9]])
    end
end

ticklabel = ["SimpleNLModels","AMPL","JuMP"]
barlabel = ["Model Creation + Solution Time (Wall)" "Solution Time (CPU)" "NLP Function Evaluation Time (CPU)"]

plt1 = groupedbar([result[1][1] result[1][2] result[1][3]], title="Luksan & Vlcek (1999)", bar_position = :dodge,
                  bar_width=0.7, xticks=(1:3, ticklabel), label=barlabel, ylabel="Time (sec)", framestyle=:box, ylim = (0,1.1*maximum(vcat(result[1]...))))
plt2 = groupedbar([result[2][1] result[2][2] result[2][3]], title="Hehn & Andrea (2011, modified)", bar_position = :dodge,
                  bar_width=0.7, xticks=(1:3, ticklabel), label=barlabel, ylabel="Time (sec)", framestyle=:box, ylim = (0,1.1*maximum(vcat(result[2]...))))
savefig(plt1,"benchmark-1.pdf")
savefig(plt2,"benchmark-2.pdf")



