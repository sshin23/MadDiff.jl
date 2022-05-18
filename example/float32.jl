# gfortran -fPIC -c -O3 -o ma27/ma27d.o ma27/ma27d.f
# gfortran -fPIC -c -O3 -o ma27/ma27s.o ma27/ma27s.f
# gfortran -fPIC -c -O3 -o common/deps.o common/deps.f
# gfortran -olibhsl.so -shared -fPIC -O3 common/deps.o ma27/ma27d.o ma27/ma27s.o -fopenmp
# $with_metis $with_blas

using MadDiff, MadNLP, MadNLPHSL, MadNLPMumps, MadNLPGPU, SparseArrays

T = Float32
N= 20;
m = MadDiffModel{T}(); 
x = [variable(m; start=mod(i,2)==1 ? -1.2 : 1.) for i=1:N]
objective(m, sum(100(x[i-1]^2-x[i])^2+(x[i-1]-1)^2 for i=2:N))
for i=1:N-2
    constraint(m, 3x[i+1]^3+2*x[i+2]-5+sin(x[i+1]-x[i+2])sin(x[i+1]+x[i+2])+4x[i+1]-x[i]exp(x[i]-x[i+1])-3 == 0)
end

import MadNLP: jac_dense!, hess_dense!


@inline function MadNLP.jac_dense!(m::MadDiffModel,x::AbstractVector,J)
    MadDiffModels.increment!(m, :neval_jac)
    MadDiffCore.jac_coord!(m.nlpcore,x,J,m.p)
    return 
end
@inline function MadNLP.hess_dense!(m::MadDiffModel,x::AbstractVector,lag::AbstractVector,z; obj_weight = 1.0)
    MadDiffModels.increment!(m, :neval_hess)
    MadDiffCore.hess_coord!(m.nlpcore,x,lag,z,m.p; obj_weight = obj_weight)
    return 
end

instantiate!(m; sparse=false) 
madnlp(
    m;
    linear_solver=MadNLPLapackGPU,
    tol = 1e-4,
    kkt_system = MadNLP.DENSE_KKT_SYSTEM,
    lapackcpu_algorithm = MadNLPLapackCPU.LU
)
# ips = MadNLP.InteriorPointSolver{
#     T,MadNLP.SparseKKTSystem{T,SparseMatrixCSC{T,Int32}}
# }(
#     m,
#     MadNLP.Options(linear_solver=MadNLPMa27, tol = 1e-4)
# )

# MadNLP.initialize!(ips.kkt)
# @time MadNLP.optimize!(ips)

