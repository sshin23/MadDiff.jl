using SimpleNLModels, Ipopt, JuMP, AmplNLWriter

# if !(@isdefined included)
    include("benchmark_include.jl")
    included = true
# end

for (N,nlm,jump) in [
    (10000,nlm_luksan_vlcek_501,jump_luksan_vlcek_501),
    (10000,nlm_ocp,jump_ocp)
]
    @time m = nlm(;N=N)
    @time instantiate!(m)
    @time optimize!(m)

    @time m = jump(;N=N,optimizer=() -> AmplNLWriter.Optimizer("ipopt"))
    @time optimize!(m)

    @time m = jump(;N=N,optimizer=Ipopt.Optimizer)
    @time optimize!(m)
end
