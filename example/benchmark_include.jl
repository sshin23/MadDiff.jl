function nlm_ocp(;N=1)
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
    
    function objective(var)
        x(i,j) = var[j+n*(i-1)]
        u(i,j) = var[j+p*(i-1)+(N+1)*n]
        return .5*sum(Q[j]*(x(i,j)-d(i,j))^2 for i=1:N for j=1:n) + .5*sum(R[j]*(u(i,j)^2) for i=1:N for j=1:p) + .5*sum(Qf[j]*(x(N+1,j)-d(N+1,j))^2 for j=1:n)
    end
    function constraint(var)
        x(i,j) = var[j+n*(i-1)]
        u(i,j) = var[j+p*(i-1)+(N+1)*n]
        return vcat(
            [-x(1,i) + x0[i] for i=1:n],
            [-x(i+1,1) + x(i,1) + (x(i,2))*dt for i=1:N],
            [-x(i+1,2) + x(i,2) + (u(i,1)*cos(x(i,7))*sin(x(i,8))*cos(x(i,9))+u(i,1)*sin(x(i,7))*sin(x(i,9)))*dt for i=1:N],
            [-x(i+1,3) + x(i,3) + (x(i,4))*dt for i=1:N], 
            [-x(i+1,4) + x(i,4) + (u(i,1)*cos(x(i,7))*sin(x(i,8))*sin(x(i,9))-u(i,1)*sin(x(i,7))*cos(x(i,9)))*dt for i=1:N], 
            [-x(i+1,5) + x(i,5) + (x(i,6))*dt for i=1:N], 
            [-x(i+1,6) + x(i,6) + (u(i,1)*cos(x(i,7))*cos(x(i,8))-9.8)*dt for i=1:N], 
            [-x(i+1,7) + x(i,7) + (b*u(i,2)*cos(x(i,7))/cos(x(i,8))+u(i,3)*sin(x(i,7))/cos(x(i,8)))*dt for i=1:N], 
            [-x(i+1,8) + x(i,8) + (-b*u(i,2)*sin(x(i,7))+u(i,3)*cos(x(i,7)))*dt for i=1:N], 
            [-x(i+1,9) + x(i,9) + (b*u(i,2)*cos(x(i,7))*tan(x(i,8))+u(i,3)*sin(x(i,7))*tan(x(i,8))+u(i,4))*dt for i=1:N]
        )
    end
    
    return IpoptProblem(objective,constraint)
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
    @NLparameter(m,d[i=1:N+1,j=1:nd]== d(i,j))
    @constraint(m,[i=1:n],x[1,i]==x0[i])
    @NLconstraint(
        m,[i=1:N],x[i+1,1] == x[i,1] + (x[i,2])*dt)
    @NLconstraint(
        m,[i=1:N], x[i+1,2] == x[i,2] + 
        (u[i,1]*cos(x[i,7])*sin(x[i,8])*cos(x[i,9])+u[i,1]*sin(x[i,7])*sin(x[i,9]))*dt)
    @NLconstraint(
        m,[i=1:N], x[i+1,3] == x[i,3] + 
        (x[i,4])*dt)
    @NLconstraint(
        m,[i=1:N], x[i+1,4] == x[i,4] + 
        (u[i,1]*cos(x[i,7])*sin(x[i,8])*sin(x[i,9])-u[i,1]*sin(x[i,7])*cos(x[i,9]))*dt)
    @NLconstraint(
        m,[i=1:N], x[i+1,5] == x[i,5] + 
        (x[i,6])*dt)
    @NLconstraint(
        m,[i=1:N], x[i+1,6] == x[i,6] + 
        (u[i,1]*cos(x[i,7])*cos(x[i,8])-9.8)*dt)
    @NLconstraint(
        m,[i=1:N], x[i+1,7] == x[i,7] + 
        (b*u[i,2]*cos(x[i,7])/cos(x[i,8])+u[i,3]*sin(x[i,7])/cos(x[i,8]))*dt)
    @NLconstraint(
        m,[i=1:N], x[i+1,8] == x[i,8] + 
        (-b*u[i,2]*sin(x[i,7])+u[i,3]*cos(x[i,7]))*dt)
    @NLconstraint(
        m,[i=1:N], x[i+1,9] == x[i,9] + 
        (b*u[i,2]*cos(x[i,7])*tan(x[i,8])+u[i,3]*sin(x[i,7])*tan(x[i,8])+u[i,4])*dt)

    @NLobjective(m,Min, .5*sum(Q[j]*(x[i,j]-d[i,j])^2 for i=1:N for j=1:n)
                 + .5*sum(R[j]*(u[i,j]^2) for i=1:N for j=1:p)
                 + .5*sum(Qf[j]*(x[N+1,j]-d[N+1,j])^2 for j=1:n))
    return m
end


function nlm_luksan_vlcek_501(;N=1)
    objective(x) = sum(100(x[i-1]^2-x[i])^2+(x[i-1]-1)^2 for i=2:N)
    constraint(x) = [
        3x[i+1]^3+2*x[i+2]-5+sin(x[i+1]-x[i+2])sin(x[i+1]+x[i+2])+4x[i+1]-x[i]exp(x[i]-x[i+1])-3
        for i=1:N-2
    ]
    x = [mod(i,2)==1 ? -1.2 : 1. for i=1:N]
    IpoptProblem(objective,constraint;x=x)
end


function jump_luksan_vlcek_501(;N=1,optimizer=Ipopt.Optimizer)
    m=Model(optimizer)
    @variable(m,x[i=1:N], start= mod(i,2)==1 ? -1.2 : 1.)
    @NLconstraint(m,[i=1:N-2], 3x[i+1]^3+2x[i+2]-5+sin(x[i+1]-x[i+2])sin(x[i+1]+x[i+2])+4x[i+1]
                -x[i]exp(x[i]-x[i+1])-3==0.)
    @NLobjective(m,Min,sum(100(x[i-1]^2-x[i])^2+(x[i-1]-1)^2 for i=2:N))
    return m
end

function nlm_eigmina(;N=1,optimizer=Ipopt.Optimizer)
    objective(x) = x[101]
    constraint(x) = [
        x[1]*x[1] + x[2]*x[2] + x[3]*x[3] + x[4]*x[4] + x[5]*x[5] + x[6]*x[6] +
        x[7]*x[7] + x[8]*x[8] + x[9]*x[9] + x[10]*x[10] + x[11]*x[11] + x[12]*x[12] +
        x[13]*x[13] + x[14]*x[14] + x[15]*x[15] + x[16]*x[16] + x[17]*x[17] + x[18]*x[18] +
        x[19]*x[19] + x[20]*x[20] + x[21]*x[21] + x[22]*x[22] + x[23]*x[23] + x[24]*x[24] +
        x[25]*x[25] + x[26]*x[26] + x[27]*x[27] + x[28]*x[28] + x[29]*x[29] + x[30]*x[30] +
        x[31]*x[31] + x[32]*x[32] + x[33]*x[33] + x[34]*x[34] + x[35]*x[35] + x[36]*x[36] +
        x[37]*x[37] + x[38]*x[38] + x[39]*x[39] + x[40]*x[40] + x[41]*x[41] + x[42]*x[42] +
        x[43]*x[43] + x[44]*x[44] + x[45]*x[45] + x[46]*x[46] + x[47]*x[47] + x[48]*x[48] +
        x[49]*x[49] + x[50]*x[50] + x[51]*x[51] + x[52]*x[52] + x[53]*x[53] + x[54]*x[54] +
        x[55]*x[55] + x[56]*x[56] + x[57]*x[57] + x[58]*x[58] + x[59]*x[59] + x[60]*x[60] +
        x[61]*x[61] + x[62]*x[62] + x[63]*x[63] + x[64]*x[64] + x[65]*x[65] + x[66]*x[66] +
        x[67]*x[67] + x[68]*x[68] + x[69]*x[69] + x[70]*x[70] + x[71]*x[71] + x[72]*x[72] +
        x[73]*x[73] + x[74]*x[74] + x[75]*x[75] + x[76]*x[76] + x[77]*x[77] + x[78]*x[78] +
        x[79]*x[79] + x[80]*x[80] + x[81]*x[81] + x[82]*x[82] + x[83]*x[83] + x[84]*x[84] +
        x[85]*x[85] + x[86]*x[86] + x[87]*x[87] + x[88]*x[88] + x[89]*x[89] + x[90]*x[90] +
        x[91]*x[91] + x[92]*x[92] + x[93]*x[93] + x[94]*x[94] + x[95]*x[95] + x[96]*x[96] +
        x[97]*x[97] + x[98]*x[98] + x[99]*x[99] + x[100]*x[100] - 1,
        x[1]*x[101] - x[1],
        x[2]*x[101] - 2*x[2],
        x[3]*x[101] - 3*x[3],
        x[4]*x[101] - 4*x[4],
        x[5]*x[101] - 5*x[5],
        x[6]*x[101] - 6*x[6],
        x[7]*x[101] - 7*x[7],
        x[8]*x[101] - 8*x[8],
        x[9]*x[101] - 9*x[9],
        x[10]*x[101] - 10*x[10],
        x[11]*x[101] - 11*x[11],
        x[12]*x[101] - 12*x[12],
        x[13]*x[101] - 13*x[13],
        x[14]*x[101] - 14*x[14],
        x[15]*x[101] - 15*x[15],
        x[16]*x[101] - 16*x[16],
        x[17]*x[101] - 17*x[17],
        x[18]*x[101] - 18*x[18],
        x[19]*x[101] - 19*x[19],
        x[20]*x[101] - 20*x[20],
        x[21]*x[101] - 21*x[21],
        x[22]*x[101] - 22*x[22],
        x[23]*x[101] - 23*x[23],
        x[24]*x[101] - 24*x[24],
        x[25]*x[101] - 25*x[25],
        x[26]*x[101] - 26*x[26],
        x[27]*x[101] - 27*x[27],
        x[28]*x[101] - 28*x[28],
        x[29]*x[101] - 29*x[29],
        x[30]*x[101] - 30*x[30],
        x[31]*x[101] - 31*x[31],
        x[32]*x[101] - 32*x[32],
        x[33]*x[101] - 33*x[33],
        x[34]*x[101] - 34*x[34],
        x[35]*x[101] - 35*x[35],
        x[36]*x[101] - 36*x[36],
        x[37]*x[101] - 37*x[37],
        x[38]*x[101] - 38*x[38],
        x[39]*x[101] - 39*x[39],
        x[40]*x[101] - 40*x[40],
        x[41]*x[101] - 41*x[41],
        x[42]*x[101] - 42*x[42],
        x[43]*x[101] - 43*x[43],
        x[44]*x[101] - 44*x[44],
        x[45]*x[101] - 45*x[45],
        x[46]*x[101] - 46*x[46],
        x[47]*x[101] - 47*x[47],
        x[48]*x[101] - 48*x[48],
        x[49]*x[101] - 49*x[49],
        x[50]*x[101] - 50*x[50],
        x[51]*x[101] - 51*x[51],
        x[52]*x[101] - 52*x[52],
        x[53]*x[101] - 53*x[53],
        x[54]*x[101] - 54*x[54],
        x[55]*x[101] - 55*x[55],
        x[56]*x[101] - 56*x[56],
        x[57]*x[101] - 57*x[57],
        x[58]*x[101] - 58*x[58],
        x[59]*x[101] - 59*x[59],
        x[60]*x[101] - 60*x[60],
        x[61]*x[101] - 61*x[61],
        x[62]*x[101] - 62*x[62],
        x[63]*x[101] - 63*x[63],
        x[64]*x[101] - 64*x[64],
        x[65]*x[101] - 65*x[65],
        x[66]*x[101] - 66*x[66],
        x[67]*x[101] - 67*x[67],
        x[68]*x[101] - 68*x[68],
        x[69]*x[101] - 69*x[69],
        x[70]*x[101] - 70*x[70],
        x[71]*x[101] - 71*x[71],
        x[72]*x[101] - 72*x[72],
        x[73]*x[101] - 73*x[73],
        x[74]*x[101] - 74*x[74],
        x[75]*x[101] - 75*x[75],
        x[76]*x[101] - 76*x[76],
        x[77]*x[101] - 77*x[77],
        x[78]*x[101] - 78*x[78],
        x[79]*x[101] - 79*x[79],
        x[80]*x[101] - 80*x[80],
        x[81]*x[101] - 81*x[81],
        x[82]*x[101] - 82*x[82],
        x[83]*x[101] - 83*x[83],
        x[84]*x[101] - 84*x[84],
        x[85]*x[101] - 85*x[85],
        x[86]*x[101] - 86*x[86],
        x[87]*x[101] - 87*x[87],
        x[88]*x[101] - 88*x[88],
        x[89]*x[101] - 89*x[89],
        x[90]*x[101] - 90*x[90],
        x[91]*x[101] - 91*x[91],
        x[92]*x[101] - 92*x[92],
        x[93]*x[101] - 93*x[93],
        x[94]*x[101] - 94*x[94],
        x[95]*x[101] - 95*x[95],
        x[96]*x[101] - 96*x[96],
        x[97]*x[101] - 97*x[97],
        x[98]*x[101] - 98*x[98],
        x[99]*x[101] - 99*x[99],
        x[100]*x[101] - 100*x[100],
    ]
    
    x = [.1 for i=1:101]
    xl= [-1. for i=1:101]
    xu= [1. for i=1:101]
    
    IpoptProblem(objective,constraint;x=x,xl=xl,xu=xu)
end

function jump_eigmina(;N=1,optimizer=Ipopt.Optimizer)
    m=Model(optimizer)
    @variable(m,-1 <= x[1:101] <= 1,start = .1)
    @constraint(m, x[1]*x[1] + x[2]*x[2] + x[3]*x[3] + x[4]*x[4] + x[5]*x[5] + x[6]*x[6] +
                x[7]*x[7] + x[8]*x[8] + x[9]*x[9] + x[10]*x[10] + x[11]*x[11] + x[12]*x[12] +
                x[13]*x[13] + x[14]*x[14] + x[15]*x[15] + x[16]*x[16] + x[17]*x[17] + x[18]*x[18] +
                x[19]*x[19] + x[20]*x[20] + x[21]*x[21] + x[22]*x[22] + x[23]*x[23] + x[24]*x[24] +
                x[25]*x[25] + x[26]*x[26] + x[27]*x[27] + x[28]*x[28] + x[29]*x[29] + x[30]*x[30] +
                x[31]*x[31] + x[32]*x[32] + x[33]*x[33] + x[34]*x[34] + x[35]*x[35] + x[36]*x[36] +
                x[37]*x[37] + x[38]*x[38] + x[39]*x[39] + x[40]*x[40] + x[41]*x[41] + x[42]*x[42] +
                x[43]*x[43] + x[44]*x[44] + x[45]*x[45] + x[46]*x[46] + x[47]*x[47] + x[48]*x[48] +
                x[49]*x[49] + x[50]*x[50] + x[51]*x[51] + x[52]*x[52] + x[53]*x[53] + x[54]*x[54] +
                x[55]*x[55] + x[56]*x[56] + x[57]*x[57] + x[58]*x[58] + x[59]*x[59] + x[60]*x[60] +
                x[61]*x[61] + x[62]*x[62] + x[63]*x[63] + x[64]*x[64] + x[65]*x[65] + x[66]*x[66] +
                x[67]*x[67] + x[68]*x[68] + x[69]*x[69] + x[70]*x[70] + x[71]*x[71] + x[72]*x[72] +
                x[73]*x[73] + x[74]*x[74] + x[75]*x[75] + x[76]*x[76] + x[77]*x[77] + x[78]*x[78] +
                x[79]*x[79] + x[80]*x[80] + x[81]*x[81] + x[82]*x[82] + x[83]*x[83] + x[84]*x[84] +
                x[85]*x[85] + x[86]*x[86] + x[87]*x[87] + x[88]*x[88] + x[89]*x[89] + x[90]*x[90] +
                x[91]*x[91] + x[92]*x[92] + x[93]*x[93] + x[94]*x[94] + x[95]*x[95] + x[96]*x[96] +
                x[97]*x[97] + x[98]*x[98] + x[99]*x[99] + x[100]*x[100] == 1)
    @constraint(m, x[1]*x[101] - x[1] == 0)
    @constraint(m, x[2]*x[101] - 2*x[2] == 0)
    @constraint(m, x[3]*x[101] - 3*x[3] == 0)
    @constraint(m, x[4]*x[101] - 4*x[4] == 0)
    @constraint(m, x[5]*x[101] - 5*x[5] == 0)
    @constraint(m, x[6]*x[101] - 6*x[6] == 0)
    @constraint(m, x[7]*x[101] - 7*x[7] == 0)
    @constraint(m, x[8]*x[101] - 8*x[8] == 0)
    @constraint(m, x[9]*x[101] - 9*x[9] == 0)
    @constraint(m, x[10]*x[101] - 10*x[10] == 0)
    @constraint(m, x[11]*x[101] - 11*x[11] == 0)
    @constraint(m, x[12]*x[101] - 12*x[12] == 0)
    @constraint(m, x[13]*x[101] - 13*x[13] == 0)
    @constraint(m, x[14]*x[101] - 14*x[14] == 0)
    @constraint(m, x[15]*x[101] - 15*x[15] == 0)
    @constraint(m, x[16]*x[101] - 16*x[16] == 0)
    @constraint(m, x[17]*x[101] - 17*x[17] == 0)
    @constraint(m, x[18]*x[101] - 18*x[18] == 0)
    @constraint(m, x[19]*x[101] - 19*x[19] == 0)
    @constraint(m, x[20]*x[101] - 20*x[20] == 0)
    @constraint(m, x[21]*x[101] - 21*x[21] == 0)
    @constraint(m, x[22]*x[101] - 22*x[22] == 0)
    @constraint(m, x[23]*x[101] - 23*x[23] == 0)
    @constraint(m, x[24]*x[101] - 24*x[24] == 0)
    @constraint(m, x[25]*x[101] - 25*x[25] == 0)
    @constraint(m, x[26]*x[101] - 26*x[26] == 0)
    @constraint(m, x[27]*x[101] - 27*x[27] == 0)
    @constraint(m, x[28]*x[101] - 28*x[28] == 0)
    @constraint(m, x[29]*x[101] - 29*x[29] == 0)
    @constraint(m, x[30]*x[101] - 30*x[30] == 0)
    @constraint(m, x[31]*x[101] - 31*x[31] == 0)
    @constraint(m, x[32]*x[101] - 32*x[32] == 0)
    @constraint(m, x[33]*x[101] - 33*x[33] == 0)
    @constraint(m, x[34]*x[101] - 34*x[34] == 0)
    @constraint(m, x[35]*x[101] - 35*x[35] == 0)
    @constraint(m, x[36]*x[101] - 36*x[36] == 0)
    @constraint(m, x[37]*x[101] - 37*x[37] == 0)
    @constraint(m, x[38]*x[101] - 38*x[38] == 0)
    @constraint(m, x[39]*x[101] - 39*x[39] == 0)
    @constraint(m, x[40]*x[101] - 40*x[40] == 0)
    @constraint(m, x[41]*x[101] - 41*x[41] == 0)
    @constraint(m, x[42]*x[101] - 42*x[42] == 0)
    @constraint(m, x[43]*x[101] - 43*x[43] == 0)
    @constraint(m, x[44]*x[101] - 44*x[44] == 0)
    @constraint(m, x[45]*x[101] - 45*x[45] == 0)
    @constraint(m, x[46]*x[101] - 46*x[46] == 0)
    @constraint(m, x[47]*x[101] - 47*x[47] == 0)
    @constraint(m, x[48]*x[101] - 48*x[48] == 0)
    @constraint(m, x[49]*x[101] - 49*x[49] == 0)
    @constraint(m, x[50]*x[101] - 50*x[50] == 0)
    @constraint(m, x[51]*x[101] - 51*x[51] == 0)
    @constraint(m, x[52]*x[101] - 52*x[52] == 0)
    @constraint(m, x[53]*x[101] - 53*x[53] == 0)
    @constraint(m, x[54]*x[101] - 54*x[54] == 0)
    @constraint(m, x[55]*x[101] - 55*x[55] == 0)
    @constraint(m, x[56]*x[101] - 56*x[56] == 0)
    @constraint(m, x[57]*x[101] - 57*x[57] == 0)
    @constraint(m, x[58]*x[101] - 58*x[58] == 0)
    @constraint(m, x[59]*x[101] - 59*x[59] == 0)
    @constraint(m, x[60]*x[101] - 60*x[60] == 0)
    @constraint(m, x[61]*x[101] - 61*x[61] == 0)
    @constraint(m, x[62]*x[101] - 62*x[62] == 0)
    @constraint(m, x[63]*x[101] - 63*x[63] == 0)
    @constraint(m, x[64]*x[101] - 64*x[64] == 0)
    @constraint(m, x[65]*x[101] - 65*x[65] == 0)
    @constraint(m, x[66]*x[101] - 66*x[66] == 0)
    @constraint(m, x[67]*x[101] - 67*x[67] == 0)
    @constraint(m, x[68]*x[101] - 68*x[68] == 0)
    @constraint(m, x[69]*x[101] - 69*x[69] == 0)
    @constraint(m, x[70]*x[101] - 70*x[70] == 0)
    @constraint(m, x[71]*x[101] - 71*x[71] == 0)
    @constraint(m, x[72]*x[101] - 72*x[72] == 0)
    @constraint(m, x[73]*x[101] - 73*x[73] == 0)
    @constraint(m, x[74]*x[101] - 74*x[74] == 0)
    @constraint(m, x[75]*x[101] - 75*x[75] == 0)
    @constraint(m, x[76]*x[101] - 76*x[76] == 0)
    @constraint(m, x[77]*x[101] - 77*x[77] == 0)
    @constraint(m, x[78]*x[101] - 78*x[78] == 0)
    @constraint(m, x[79]*x[101] - 79*x[79] == 0)
    @constraint(m, x[80]*x[101] - 80*x[80] == 0)
    @constraint(m, x[81]*x[101] - 81*x[81] == 0)
    @constraint(m, x[82]*x[101] - 82*x[82] == 0)
    @constraint(m, x[83]*x[101] - 83*x[83] == 0)
    @constraint(m, x[84]*x[101] - 84*x[84] == 0)
    @constraint(m, x[85]*x[101] - 85*x[85] == 0)
    @constraint(m, x[86]*x[101] - 86*x[86] == 0)
    @constraint(m, x[87]*x[101] - 87*x[87] == 0)
    @constraint(m, x[88]*x[101] - 88*x[88] == 0)
    @constraint(m, x[89]*x[101] - 89*x[89] == 0)
    @constraint(m, x[90]*x[101] - 90*x[90] == 0)
    @constraint(m, x[91]*x[101] - 91*x[91] == 0)
    @constraint(m, x[92]*x[101] - 92*x[92] == 0)
    @constraint(m, x[93]*x[101] - 93*x[93] == 0)
    @constraint(m, x[94]*x[101] - 94*x[94] == 0)
    @constraint(m, x[95]*x[101] - 95*x[95] == 0)
    @constraint(m, x[96]*x[101] - 96*x[96] == 0)
    @constraint(m, x[97]*x[101] - 97*x[97] == 0)
    @constraint(m, x[98]*x[101] - 98*x[98] == 0)
    @constraint(m, x[99]*x[101] - 99*x[99] == 0)
    @constraint(m, x[100]*x[101] - 100*x[100] == 0)
    @objective(m, Min, x[101])
    return m
end
