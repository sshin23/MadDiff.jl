import NonlinearModels: IpoptProblem
using NonlinearModels, MadNLP, Ipopt, JuMP, AmplNLWriter

include("benchmark_include.jl")

for (N,nlm,jump) in [
    (10000,nlm_luksan_vlcek_501,jump_luksan_vlcek_501),
    (10000,nlm_ocp,jump_ocp),
    (1,nlm_eigmina,jump_eigmina),
]
    @time iprob = nlm(;N=N)
    @time Ipopt.solveProblem(iprob)

    @time m = jump(;N=N,optimizer=() -> AmplNLWriter.Optimizer("ipopt"))
    @time optimize!(m)

    @time m = jump(;N=N,optimizer=Ipopt.Optimizer)
    @time optimize!(m)
end
