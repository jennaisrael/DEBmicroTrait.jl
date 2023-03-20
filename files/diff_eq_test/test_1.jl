#Call necessary packages from Julia
using DifferentialEquations
using Plots
using LaTeXStrings
function plume!(du,u,p,t) 
    # let u be the vector of variables [Q,M,B], p be a vector of constants, and t be the time vector
    α, Nsqr = p # the constants alpha and N^2 as defined in our system
    du[1] = 2.0*α.*pi.^(0.5).*(u[2].^(0.5)) #dQ/dz in terms of M
    du[2] = u[3].*u[1]./u[2]            #dM/dz in terms of B, Q, and M
    du[3] = -1*Nsqr.*u[1]                  #dB/dz in terms of Q
end

l = @layout [a; b; c] #initialize subplot layout


# Case 1: N^2 = 0, M0=0
p = (0.1, 0.0)# parameter collection alpha and N2
u0 = [0.00001; 0.00001; 2.0] # system is a bit stiff at the beginning, using these u0 values > eps seems to work okay
tspan = (0.0,10.0) # zmax=10 m
prob = ODEProblem(plume!,u0,tspan,p) #sets up the system, intial values, range of z values to calculate over, and constants to be passed to the system
sol = solve(prob, Tsit5()) #runs the solver using a fourth-order, five-stage explicit Runge-Kutta algorithm
#sol output has properties t, or the vector of z values, and u, the numerically calculated values of Q, M, and B at each z coordinate in t
w=sol[2,:]./sol[1,:] #w bar calculated as M/Q

#Unstratified ambient analytical solution
C=(6*0.1/5)*(9*0.1/10)^(1/3)*pi^(2/3)
Q=C*2.0^(1/3)*sol.t.^(5/3)
M=(9*0.1/10)^(2/3)*pi^(1/3)*2.0^(2/3)*sol.t.^(4/3)
B=2.0*ones(size(sol.t))

#Plot case 1
p1=plot([sol[1,:] sol[2,:] sol[3,:] w Q M B M./Q], sol.t, title= L"Case 1", label=[L"Q [m^3/s]" L"M [m^{4}/s^{2}]" L"B [m^{4}/s^{3}]" L"\bar{w} [m/s]" L"Q_{analytical}" L"M_{analytical}" L"B_{analytical}" L"\bar{w}_{analytical}" ], linewidth=3, ls=[:solid :solid :solid :solid :dash :dash :dash :dash])
xlims!((0,15)) #limit the x values displayed bc wbar looks weird at the beginning
ylabel!(L"z [m]")


#########################################################
#Case 2: N^2 = 0, M0=1 m^4/s^2
p = (0.1, 0.0) # parameter collection alpha and N2
u0 = [0.0;1.0;2.0] #Q0, M0, B0
tspan = (0.0,10.0) # zmax=10 m
prob = ODEProblem(plume!,u0,tspan,p) #sets up the system, intial values, range of z values to calculate over, and constants to be passed to the system
sol = solve(prob, Tsit5()) #runs the solver using a fourth-order, five-stage explicit Runge-Kutta algorithm
w=sol[2,:]./sol[1,:] #w bar calculated as M/Q

#Need to re-calculate the analytical solution for each case because the time (z) vector is allocated differently each time
C=(6*0.1/5)*(9*0.1/10)^(1/3)*pi^(2/3)
Q=C*2.0^(1/3)*sol.t.^(5/3)
M=(9*0.1/10)^(2/3)*pi^(1/3)*2.0^(2/3)*sol.t.^(4/3)
B=2.0*ones(size(sol.t))

#Plot Case 2
p2 = plot( [sol[1,:] sol[2,:] sol[3,:] w Q M B M./Q], sol.t, title= L"Case 2", label=[L"Q [m^3/s]" L"M [m^{4}/s^{2}]" L"B [m^{4}/s^{3}]" L"\bar{w} [m/s]" L"Q_{analytical}" L"M_{analytical}" L"B_{analytical}" L"\bar{w}_{analytical}" ], linewidth=3, ls=[:solid :solid :solid :solid :dash :dash :dash :dash])
xlims!((0,15))
ylabel!(L"z [m]")

############################################################
#Case 3: N^2 = 0.075s^-2, M=0
p = (0.1, 0.075) # parameter collection alpha and N2
u0 = [1.0*10^(-5);1.0*10^(-6);2.0] #Q0, M0, B0
tspan = (0.0,10.0) # zmax=10 m
prob = ODEProblem(plume!,u0,tspan,p)
# Try a differnt solver, reccommended for stiff problems, second order A-B-L-S-stable one-step ESDIRK method
sol = solve(prob, TRBDF2())

#Need to re-calculate the analytical solution for each case because the time (z) vector is allocated differently each time
w=sol[2,:]./sol[1,:]
C=(6*0.1/5)*(9*0.1/10)^(1/3)*pi^(2/3)
Q=C*2.0^(1/3)*sol.t.^(5/3)
M=(9*0.1/10)^(2/3)*pi^(1/3)*2.0^(2/3)*sol.t.^(4/3)
B=2.0*ones(size(sol.t))

#Plot case 3
p3 = plot( [sol[1,:] sol[2,:] sol[3,:] w Q M B M./Q], sol.t, title= L"Case 3", label=[L"Q [m^3/s]" L"M [m^{4}/s^{2}]" L"B [m^{4}/s^{3}]" L"\bar{w} [m/s]" L"Q_{analytical}" L"M_{analytical}" L"B_{analytical}" L"\bar{w}_{analytical}" ], linewidth=3, ls=[:solid :solid :solid :solid :dash :dash :dash :dash])
xlims!((0,15))
ylabel!(L"z [m]")

plot(p1, p2, p3, layout = l)


# ## Example from documentation
# # using DifferentialEquations
# # function lorenz!(du,u,p,t)
# #     σ,ρ,β = p
# #     du[1] = σ*(u[2]-u[1])
# #     du[2] = u[1]*(ρ-u[3]) - u[2]
# #     du[3] = u[1]*u[2] - β*u[3]
# # end
# # u0 = [1.0,0.0,0.0]
# # p = (10,28,8/3)
# # tspan = (0.0,100.0)
# # prob = ODEProblem(lorenz!,u0,tspan,p)
# # sol = solve(prob)
# # plot(sol)