function nlm_ocp(;N=1,optimizer=Ipopt.Optimizer,opt...)
    q = 1
    b = 1
    n = 9
    p = 4
    nd= 9
    x0 = [0,0,0,0,0,0,0,0,0]
    dt = .01
    Q = [1,0,1,0,1,0,1,1,1]
    Qf= [1,0,1,0,1,0,1,1,1]/dt
    R = ones(4)/10
    
    d(i,j) = (j==1 ? 1*sin(2*pi/N*i) : 0) + (j==3 ? 2*sin(4*pi/N*i) : 0) + (j==5 ? 2*i/N : 0)

    m = SimpleNLModels.Model(optimizer;opt...)
    
    x=[variable(m,start = 0) for i=1:N+1,j=1:n]
    u=[variable(m,start = 0) for i=1:N,j=1:p]
    [constraint(m,x[1,i]-x0[i]) for i=1:n]
    for i=1:N
        constraint(m, -x[i+1,1] + x[i,1] + (x[i,2])*dt)
        constraint(m, -x[i+1,2] + x[i,2] + (u[i,1]*cos(x[i,7])*sin(x[i,8])*cos(x[i,9])+u[i,1]*sin(x[i,7])*sin(x[i,9]))*dt)
        constraint(m, -x[i+1,3] + x[i,3] + (x[i,4])*dt)
        constraint(m, -x[i+1,4] + x[i,4] + (u[i,1]*cos(x[i,7])*sin(x[i,8])*sin(x[i,9])-u[i,1]*sin(x[i,7])*cos(x[i,9]))*dt)
        constraint(m, -x[i+1,5] + x[i,5] + (x[i,6])*dt)
        constraint(m, -x[i+1,6] + x[i,6] + (u[i,1]*cos(x[i,7])*cos(x[i,8])-9.8)*dt)
        constraint(m, -x[i+1,7] + x[i,7] + (b*u[i,2]*cos(x[i,7])/cos(x[i,8])+u[i,3]*sin(x[i,7])/cos(x[i,8]))*dt)
        constraint(m, -x[i+1,8] + x[i,8] + (-b*u[i,2]*sin(x[i,7])+u[i,3]*cos(x[i,7]))*dt)
        constraint(m, -x[i+1,9] + x[i,9] + (b*u[i,2]*cos(x[i,7])*tan(x[i,8])+u[i,3]*sin(x[i,7])*tan(x[i,8])+u[i,4])*dt)
    end
    for i=1:N
        for j=1:n
            objective(m,.5*Q[j]*(x[i,j]-d(i,j))^2 )
        end
        for j=1:p
            objective(m,.5*R[j]*(u[i,j]^2))
        end
    end
    for j=1:n
        objective(m,.5*Qf[j]*(x[N+1,j]-d(N+1,j))^2)
    end
        
    return m
end

function jump_ocp(;N=1,optimizer=Ipopt.Optimizer)
    q = 1
    b = 1
    n = 9
    p = 4
    nd= 9
    x0 = [0,0,0,0,0,0,0,0,0]
    d = (i,j)->
        (j==1 ? 1*sin(2*pi/N*i) : 0) +
        (j==3 ? 2*sin(4*pi/N*i) : 0) +
        (j==5 ? 2*i/N : 0)
    dt = .01

    Q = [1,0,1,0,1,0,1,1,1]
    Qf= [1,0,1,0,1,0,1,1,1]/dt
    R = ones(4)/10
    
    m = Model(optimizer)
    @variable(m,x[1:N+1,1:n],start = 0)
    @variable(m,u[1:N,1:p],start = 0)
    @constraint(m,[i=1:n],x[1,i]==x0[i])
    @constraint(m,[i=1:N],x[i+1,1] == x[i,1] + (x[i,2])*dt)
    @NLconstraint(m,[i=1:N], x[i+1,2] == x[i,2] + (u[i,1]*cos(x[i,7])*sin(x[i,8])*cos(x[i,9])+u[i,1]*sin(x[i,7])*sin(x[i,9]))*dt)
    @constraint(m,[i=1:N], x[i+1,3] == x[i,3] + (x[i,4])*dt)
    @NLconstraint(m,[i=1:N], x[i+1,4] == x[i,4] + (u[i,1]*cos(x[i,7])*sin(x[i,8])*sin(x[i,9])-u[i,1]*sin(x[i,7])*cos(x[i,9]))*dt)
    @constraint(m,[i=1:N], x[i+1,5] == x[i,5] + (x[i,6])*dt)
    @NLconstraint(m,[i=1:N], x[i+1,6] == x[i,6] + (u[i,1]*cos(x[i,7])*cos(x[i,8])-9.8)*dt)
    @NLconstraint(m,[i=1:N], x[i+1,7] == x[i,7] + (b*u[i,2]*cos(x[i,7])/cos(x[i,8])+u[i,3]*sin(x[i,7])/cos(x[i,8]))*dt)
    @NLconstraint(m,[i=1:N], x[i+1,8] == x[i,8] + (-b*u[i,2]*sin(x[i,7])+u[i,3]*cos(x[i,7]))*dt)
    @NLconstraint(m,[i=1:N], x[i+1,9] == x[i,9] + (b*u[i,2]*cos(x[i,7])*tan(x[i,8])+u[i,3]*sin(x[i,7])*tan(x[i,8])+u[i,4])*dt)
    @objective(m,Min, .5*sum(Q[j]*(x[i,j]-d(i,j))^2 for i=1:N for j=1:n) + .5*sum(R[j]*(u[i,j]^2) for i=1:N for j=1:p)
                 + .5*sum(Qf[j]*(x[N+1,j]-d(N+1,j))^2 for j=1:n))
    return m
end


function nlm_luksan_vlcek_501(;N=1,optimizer=Ipopt.Optimizer,opt...)
    m = SimpleNLModels.Model(optimizer;opt...)

    x = [variable(m;start=mod(i,2)==1 ? -1.2 : 1.) for i=1:N];

    for i=2:N
        objective(m,100(x[i-1]^2-x[i])^2+(x[i-1]-1)^2)
    end

    for i=1:N-2
        constraint(m,3x[i+1]^3+2*x[i+2]-5+sin(x[i+1]-x[i+2])sin(x[i+1]+x[i+2])+4x[i+1]-x[i]exp(x[i]-x[i+1])-3)
    end
    
    return m
end

function jump_luksan_vlcek_501(;N=1,optimizer=Ipopt.Optimizer)
    m=Model(optimizer)
    @variable(m,x[i=1:N], start= mod(i,2)==1 ? -1.2 : 1.)
    @NLconstraint(m,[i=1:N-2], 3x[i+1]^3+2x[i+2]-5+sin(x[i+1]-x[i+2])sin(x[i+1]+x[i+2])+4x[i+1]-x[i]exp(x[i]-x[i+1])-3==0.)
    @NLobjective(m,Min,sum(100(x[i-1]^2-x[i])^2+(x[i-1]-1)^2 for i=2:N))
    return m
end

function parse_ipopt_output(filename)
    output = read(filename,String);
    i = match(r"Total CPU secs in IPOPT",output).offset
    j = match(r"Total CPU secs in NLP function evaluations",output).offset;
    k = match(r"EXIT:",output).offset

    cpu_time_ipopt = parse(Float64,output[i+54:j-3])
    cpu_time_nlp_function = parse(Float64,output[j+54:k-3])
    cpu_time_total = cpu_time_ipopt + cpu_time_nlp_function

    return cpu_time_total, cpu_time_nlp_function
end
