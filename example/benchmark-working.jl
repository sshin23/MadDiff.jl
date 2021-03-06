using CasADi

function nlm_luksan_vlcek_501(;N=1,optimizer=SimpleNLModels.IpoptOptimizer,opt...)

    m = casadi.Opti()
    
    x = [opti._variable() for i=1:N];
    
    for i=1:N
        opti.set_initial(x[i],mod(i,2)==1 ? -1.2 : 1.)
    end
    
    opti.minimize(sum(100(x[i-1]^2-x[i])^2+(x[i-1]-1)^2 for i=2:N))
    for i=1:N-2
        opti._subject_to(3x[i+1]^3+2*x[i+2]-5+sin(x[i+1]-x[i+2])sin(x[i+1]+x[i+2])+4x[i+1]-x[i]*2.718281828459045^(x[i]-x[i+1])-3 == 0)
    end

    
    
    return m
end
