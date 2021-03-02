get_num_variables(ex::Expression) = maximum(i for (i,d) in deriv(ex))
get_num_variables(exs::Expression...) = maximum(get_num_variables(ex) for ex in exs)

function fill!(f,pairs,x)
    for (i,d) in pairs
        f[i] += d(x)
    end
end

function fill!(J,pairs,x,offset)
    for (j,d) in pairs
        J[offset += 1] = d(x)
    end
    return offset
end

get_pairss(dict) =  [
    [(i,d) for (i,d) in dict if typeof(d)==type]
    for type in union(typeof(d) for (i,d) in dict)
]
get_con_pairss(cons) = [
    [(i,func(cons[i])) for i=1:length(cons) if typeof(func(cons[i]))==type]
    for type in union(typeof(func(con)) for con in cons)
]

function eval_grad!(obj)
    pairss = get_pairss(deriv(obj))
    function (f,x)
        f.=0
        for pairs in pairss
            fill!(f,pairs,x)
        end
    end
end

function eval_con!(cons)
    pairss = get_con_pairss(cons)
    function (c,x)
        c.= 0
        for pairs in pairss
            fill!(c,pairs,x)
        end
    end
end

function get_jac_dict(cons)
    jac_dict = Dict{Tuple{Int,Int},Function}()
    for i=1:length(cons)
        for (j,d) in deriv(cons[i])
            if haskey(jac_dict,(i,j))
                dold = jac_dict[i,j]
                jac_dict[i,j] = x->dold(x)+d(x)
            else
                jac_dict[i,j] = x->d(x)
            end
        end
    end
    return jac_dict
end

function eval_jac!(cons)
    jac_dict = get_jac_dict(cons)
    pairss = get_pairss(jac_dict)
    
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

function get_hess_dict(obj_2nd,con_2nds)
    hess_dict = Dict{Tuple{Int,Int},Function}()

    for (i,dobj) in obj_2nd
        for (j,d) in dobj
            if haskey(hess_dict,(i,j))
                dold = hess_dict[i,j]
                hess_dict[i,j] = (x,sig,lag)->dold(x,sig,lag)+d(x)*sig
            else
                i>=j && (hess_dict[i,j] = (x,sig,lag)->d(x)*sig)
            end
        end
    end

    for k=1:length(con_2nds)
        for (i,dcon) in con_2nds[k]
            for (j,d) in dcon
                if haskey(hess_dict,(i,j))
                    dold = hess_dict[i,j]
                    hess_dict[i,j] = (x,sig,lag)->dold(x,sig,lag)+d(x)*lag[k]
                else
                    i>=j && (hess_dict[i,j] = (x,sig,lag)->d(x)*lag[k])
                end
            end
        end
    end
    return hess_dict
end

get_obj_2nd(obj) = [(i,deriv(d(Source()))) for (i,d) in deriv(obj)]
get_con_2nds(cons) = [[(i,deriv(d(Source()))) for (i,d) in deriv(con)] for con in cons]

function eval_hess!(obj,cons)
    obj_2nd = get_obj_2nd(obj)
    con_2nds = get_con_2nds(cons)
    
    hess_dict = get_hess_dict(obj_2nd,con_2nds)
    pairss = get_pairss(hess_dict)
    
    function hess!(H,x,lag,sig)
        offset = 0
        for pairs in pairss
            offset = fill_hessian!(H,pairs,x,sig,lag,offset)
        end
    end

    function hess_sparsity!(I,J)
        offset = 0
        for pairs in pairss
            offset = sparsity!(I,J,pairs,offset)
        end
    end

    nnz_hess = sum(length(pairs) for pairs in pairss; init = 0) 

    (hess!,hess_sparsity!,nnz_hess)
end

function fill_hessian!(H,pairs,x,sig,lag,offset)
    for ((i,j),d) in pairs
        H[offset += 1] = d(x,sig,lag)
    end
    return offset
end

function sparsity!(I,J,pairs,offset)
    for ((i,j),d) in pairs
        offset += 1
        I[offset] = i
        J[offset] = j
    end
    return offset
end


function get_nlp_functions(obj,cons)

    _obj = func(obj)
    _grad! = eval_grad!(obj)
    _con! = eval_con!(cons)
    _jac!,_jac_sparsity!,nnz_jac = eval_jac!(cons)
    _hess!,_hess_sparsity!,nnz_hess = eval_hess!(obj,cons)
    
    (_obj,_grad!,_con!,_jac!,_jac_sparsity!,nnz_jac,_hess!,_hess_sparsity!,nnz_hess)
end
