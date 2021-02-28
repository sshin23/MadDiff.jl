get_num_variables(ex::Expression) = maximum(i for (i,d) in der(ex))
get_num_variables(exs::Expression...) = maximum(get_num_variables(ex) for ex in exs)

function fill!(f,dobj,x)
    for (i,d) in dobj
        f[i] += d(x)
    end
end

function fill!(J,dcon,x,offset)
    for (j,d) in dcon
        J[offset += 1] = d(x)
    end
    return offset
end

get_derivative_pairss(dict) =[
    [(i,d) for (i,d) in dict if typeof(d)==type]
    for type in union(typeof(d) for (i,d) in dict)
]

function eval!(pairss)
    function (f,x)
        f.=0
        for pairs in pairss
            fill!(f,pairs,x)
        end
    end
end
function eval_jac!(pairss)
    function jac!(J,x)
        offset = 0
        for pairs in pairss
            offset = fill!(J,pairs,x,offset)
        end
    end
    function jac_sparsity!(I,J)
        offset = 0
        for pairs in pairss
            offset = sparsity!(I,J,pairs,offset)
        end
    end
    nnz_jac = sum(length(pairs) for pairs in pairss; init = 0)
    
    (jac!,jac_sparsity!,nnz_jac)
end

function eval_hess!(pairss)
    function hess!(H,x,lag,sig)
        offset = 0
        for obj_pairs in pairss
            offset = fill_hessian!(H,obj_pairs,x,sig,lag,offset)
        end
    end

    function hess_sparsity!(I,J)
        offset = 0
        for obj_pairs in pairss
            offset = sparsity!(I,J,obj_pairs,offset)
        end
    end

    nnz_hess = sum(length(pairs) for pairs in pairss; init = 0) 

    (hess!,hess_sparsity!,nnz_hess)
end

function fill_hessian!(H,dict,x,sig,lag,offset)
    for ((i,j),d) in dict
        H[offset += 1] = d(x,sig,lag)
    end
    return offset
end

function sparsity!(I,J,dict,offset)
    for ((i,j),d) in dict
        offset += 1
        I[offset] = i
        J[offset] = j
    end
    return offset
end


function get_nlp_functions(obj,cons,n,m)
    obj_2nd = [i=>der(d(Expression())) for (i,d) in der(obj)]
    con_2nds = [[i=>der(d(Expression())) for (i,d) in der(con)] for con in cons]
    
    merged = Dict{Tuple{Int,Int},Function}()

    for (i,dobj) in obj_2nd
        for (j,d) in dobj
            if haskey(merged,(i,j))
                dold = merged[i,j]
                merged[i,j] = (x,sig,lag)->dold(x,sig,lag)+d(x)*sig
            else
                i>=j && (merged[i,j] = (x,sig,lag)->d(x)*sig)
            end
        end
    end

    for k=1:m 
        for (i,dcon) in con_2nds[k]
            for (j,d) in dcon
                if haskey(merged,(i,j))
                    dold = merged[i,j]
                    merged[i,j] = (x,sig,lag)->dold(x,sig,lag)+d(x)*lag[k]
                else
                    i>=j && (merged[i,j] = (x,sig,lag)->d(x)*lag[k])
                end
            end
        end
    end
    
    
    _obj = obj.fun
    _grad! = eval!(get_derivative_pairss(der(obj)))
    _con! = eval!(get_derivative_pairss([i=>fun(cons[i]) for i=1:length(cons)]))

    conjac = Dict{Tuple{Int,Int},Function}()
    for i=1:length(cons)
        for (j,d) in der(cons[i])
            if haskey(conjac,(i,j))
                dold = conjac[i,j]
                conjac[i,j] = x->dold(x)+d(x)
            else
                conjac[i,j] = x->d(x)
            end
        end
    end
    _jac!,_jac_sparsity!,nnz_jac = eval_jac!(get_derivative_pairss(conjac))
    _hess!,_hess_sparsity!,nnz_hess = eval_hess!(get_derivative_pairss(merged))
    
    (_obj,_grad!,_con!,_jac!,_jac_sparsity!,nnz_jac,_hess!,_hess_sparsity!,nnz_hess)
end
