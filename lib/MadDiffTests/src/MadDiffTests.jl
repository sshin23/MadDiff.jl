module MadDiffTests

using Test, MadDiffCore, ForwardDiff

const atol = 1e-6
const rtol = 1e-6
function compare(p1,p2)
    for i=1:length(p1)
        isequal(p1[i], NaN) && isequal(p2[i], NaN) && continue
        d = abs(p1[i] .- p2[i])
        
        d > atol && d / max(abs(p1[i]),abs(p2[i])) > rtol && return false
    end
    return true
end

function compare(f1,f2,xs,ps)
    for i=1:length(xs)
        compare(f1(xs[i],ps[i]),f2(xs[i],ps[i])) || return false
    end
    return true
end

function compare(f1,f2,y1,y2,xs,ps)
    for i=1:length(xs)
        f1(y1,xs[i],ps[i])
        f2(y2,xs[i],ps[i])
        compare(y1,y2) || return false
    end
    return true
end


function _force_lower_tri(z)
    for i=1:size(z,1)
        for j=1:size(z,2)
            i<j && (z[i,j]=0)
        end
    end
end

_get_gorig!(forig) = (y,x,p) -> (y.=0;ForwardDiff.gradient!(y,x->forig(x,p),x))
_get_sgorig!(forig,gsparsity) = (y,x,p) ->  (y.=0;y[1:length(gsparsity)] .= ForwardDiff.gradient(x->forig(x,p),x)[gsparsity])
_get_horig!(forig) = (z,x,p) -> (z.=0;ForwardDiff.hessian!(z,x->forig(x,p),x); _force_lower_tri(z))
_get_shorig!(forig,hsparsity) = (z,x,p) ->  (z.=0;hess=ForwardDiff.hessian(x->forig(x,p),x); z[1:length(hsparsity)] .= [hess[i,j] for (i,j) in hsparsity])
_get_jorig!(forig) = (z,x,p) -> (z.=0;ForwardDiff.jacobian!(z,(y,x)->forig(y,x,p),zeros(size(z,1)),x))
_get_sjorig!(forig,jsparsity) = (y,x,p) ->  (y.=0;n=length(y);jac=ForwardDiff.jacobian((y,x)->forig(y,x,p),zeros(n),x); y[1:length(jsparsity)] .= [jac[i,j] for (i,j) in jsparsity])

function test_function(forig,name,xs,ps,y1,y2,z1,z2)
    ff = forig(MadDiffCore.Variable(),MadDiffCore.Parameter())
    f  = MadDiffCore.Evaluator(ff)
    
    g! = MadDiffCore.GradientEvaluator(ff)
    sg! = MadDiffCore.SparseGradientEvaluator(ff)
    
    h! = MadDiffCore.HessianEvaluator(ff)
    sh! = MadDiffCore.SparseHessianEvaluator(ff)
    
    gorig!  = _get_gorig!(forig)
    sgorig! = _get_sgorig!(forig,sg!.sparsity)
    horig!  = _get_horig!(forig)
    shorig!  = _get_shorig!(forig,sh!.sparsity)

    @testset "$name" begin
        @test compare(forig,f,xs,ps)
        @test compare(gorig!,g!,y1,y2,xs,ps)
        @test compare(sgorig!,sg!,y1,y2,xs,ps)
        @test compare(horig!,h!,z1,z2,xs,ps) || (println(z1);println(z2);false)
        @test compare(shorig!,sh!,z1,z2,xs,ps)
    end
end

function test_field(forig,name,xs,ps,y1,y2,z1,z2)
    ff = MadDiffCore.Field()
    forig(ff,MadDiffCore.Variable(),MadDiffCore.Parameter())
    
    f  = MadDiffCore.FieldEvaluator(ff)
    
    j! = MadDiffCore.JacobianEvaluator(ff)
    sj! = MadDiffCore.SparseJacobianEvaluator(ff)
    
    jorig!  = _get_jorig!(forig)
    sjorig! = _get_sjorig!(forig,sj!.sparsity)

    @testset "$name" begin
        @test compare(forig,f,y1,y2,xs,ps)
        @test compare(jorig!,j!,z1,z2,xs,ps)
        @test compare(sjorig!,sj!,y1,y2,xs,ps)
    end
end

end # module
