using DEBmicroTrait
using CSV, DataFrames, Statistics
using JLD
using DifferentialEquations
using Plots
using LaTeXStrings

dir                     = "DEBSCRIPTS" in keys(ENV) ? ENV["DEBSCRIPTS"] : pwd()

# Load media composition: Formula, Molecular weight, Medium concentration
df_metabolites = CSV.read(joinpath(dir, "files/output2/LB_medium_unique.csv"), DataFrame, missingstring="N/A")
df_isolates             = CSV.read(joinpath(dir, "files/input/isolates2traits.csv"), DataFrame, missingstring="N/A")


# Load isolate parameterization, note: assimilation parameterization depends on media composition
assimilation            = load(joinpath(dir, "files/output//isolates_assimilation_LB.jld")) # run isolates_assimilation_1_10_R2A.jl first
enzymes                 = load(joinpath(dir, "files/output/isolates_enzymes.jld"))
maintenance             = load(joinpath(dir, "files/output/isolates_maintenance.jld"))
protein_synthesis       = load(joinpath(dir, "files/output/isolates_protein_synthesis.jld"))
turnover                = load(joinpath(dir, "files/output/isolates_turnover.jld"))
initb                   = load(joinpath(dir, "files/output/isolates_batch_init.jld"))


id_isolate = 7 #HA54 is 7, HB15 is 30
n_isolates = length(id_isolate)
n_monomers = 22 #22 unique, previously 38

l = @layout [a;b] #initialize subplot layout

#try just plotting HB15 and HA54

#plot the initial concentrations in LB vs the KD values in the assimilation dictionary
KD=get(assimilation, "KD",1) # 22 metabolites by 39 isolates
x=LinRange(1,22,22)
label=df_metabolites.Compound #label x ticks
p1=plot(x, KD[:,[7,30]], label=["HB15" "HA54"])#label=df_isolates.Abbreviation)
plot!(x,df_metabolites.Concentration_sum,label="LB C_O",linestyle=:dash,linewidth=3,lc=:black)
plot!(xticks = (x, label), xrotation=90)

ylabel!("[Substrate] (mM)")
#xlabel!("Substrate")
title!("LB")

#repeat for R2A

df_R2A = CSV.read(joinpath(dir, "files/output2/1_10_R2A_medium_unique.csv"), DataFrame, missingstring="N/A")
assimilation_R2A            = load(joinpath(dir, "files/output//isolates_assimilation_10_R2A.jld")) 

KD2=get(assimilation_R2A, "KD",1)
x2=LinRange(1,23,23)

p2=plot(x2,KD2[:,[7,30]] , label=["HB15" "HA54"] )#label=df_isolates.Abbreviation)
plot!(x2, df_R2A.Concentration_sum, label="1/10 R2A C_0",linestyle=:dash,linewidth=3,lc=:black)
plot!(xticks = (x2, df_R2A.Compound), xrotation=90)

ylabel!("[Substrate] (mM)")
#xlabel!("Substrate")
title!("1/10 R2A")

p3= plot(p1, p2, layout = l, size=(1200,800),left_margin=[20mm 0mm],bottom_margin=[20mm 0mm])
display(p3)


#make another plot with the KDs for LB and R2A on one plot, R2A has one more monomer 
#than LB (glucose)


#concatenate the HB15 and HA54 KD vector to the metabolites vector
df_metabolites[:,:HA54_KD] = KD[:,7]
df_metabolites[:,:HB15_KD] = KD[:,30]
df_R2A[:,:HA54_KD] = KD2[:,7]
df_R2A[:,:HB15_KD] = KD2[:,30]
#dataframes do not have substrates in the same order, need to order first
df_metabolites_sorted=sort!(df_metabolites,[order(:Compound)])
df_R2A_sorted = sort!(df_R2A,[order(:Compound)])

#What is the index of glucose in df_R2A_sorted? don't plot a KD value for LB at that index
df_R2A_sorted[:,:RowNum]=1:length(df_R2A_sorted.Compound)
glucose_row=filter(:Compound => ==("glucose"), df_R2A_sorted)
gind=glucose_row.RowNum[1]

p4=plot(x2,df_R2A_sorted.HA54_KD , label="HA54 in 1/10 R2A", linestyle=:dash,linewidth=3)
plot!(x2,df_R2A_sorted.HB15_KD , label="HB15 in 1/10 R2A", linestyle=:dash,linewidth=3)

plot!([x2[1:gind-1]; x2[gind+1:end]],df_metabolites_sorted.HA54_KD,label="HA54 in LB", linestyle=:dot,linewidth=3)
plot!([x2[1:gind-1]; x2[gind+1:end]],df_metabolites_sorted.HB15_KD,label="HB15 in LB", linestyle=:dot,linewidth=3)

#add some KD values calculated for one strain and one substrate
assimilationHA54gly            = load(joinpath(dir, "files/output//isolates_assimilation_HA54_glycine.jld"))
assimilationHB15ser            = load(joinpath(dir, "files/output//isolates_assimilation_HB15_serine.jld"))

KD54gly=get(assimilationHA54gly, "KD",1)
KD15ser=get(assimilationHB15ser, "KD",1)

#indices of glycine and serine
glycine_row=filter(:Compound => ==("glycine") , df_R2A_sorted)
serine_row=filter(:Compound => ==("serine") , df_R2A_sorted)
glyind=glycine_row.RowNum[1]
serind=serine_row.RowNum[1]
plot!((glyind, KD54gly[1]), seriestype=:scatter, label="HA54 in glycine")
plot!((serind, KD15ser[1]), seriestype=:scatter, label="HB15 in serine")

ylabel!("[Substrate] (mM)")
plot!(xticks = (x2, df_R2A_sorted.Compound), xrotation=90)
margin=20mm

l2= @layout[a]
p5= plot(p4, layout =l2, size=(1200,800),left_margin=[20mm 0mm],bottom_margin=[20mm 0mm])
display(p5)


# Determine if the substrate consumption profiles when run together in LB as when each is passed to the model individually. 
# Compare affinity paramters for LB run with individually run extimates.

# p                 = DEBmicroTrait.init_mixed_medium(id_isolate, n_monomers, assimilation, enzymes, maintenance, protein_synthesis, turnover)
# n_polymers        = p.setup_pars.n_polymers
# n_monomers        = p.setup_pars.n_monomers
# n_microbes        = p.setup_pars.n_microbes

#  #initial conditions for substrate, reserve, and structural biomass (should have same )
# u0                                                                         = zeros(p.setup_pars.dim) #47, 43 monomers+ 0 polymers+ 1 strain + reserve biomass+ structure biomass+ enzyme concentration +total respiration
# u0[1+n_polymers+n_monomers:n_polymers+n_monomers+n_microbes]              .= 0.9*initb["Bio0"][id_isolate] #90% of initial biomass is reserve
# u0[1+n_polymers+n_monomers+n_microbes:n_polymers+n_monomers+2*n_microbes] .= 0.1*initb["Bio0"][id_isolate] #10% structure
# u0[1+n_polymers:n_polymers+n_monomers]                                    .= df_metabolites.Concentration_sum # initalizes the concentrations of monomers 



# tspan             = (0.0,192.0) #time in hours
# prob              = ODEProblem(DEBmicroTrait.batch_model!,u0,tspan,p)
# sol               = solve(prob, alg_hints=[:stiff]) #first dimension is variable

# n_t=length(sol.t)

# D_tseries =  sol[1+n_polymers:n_polymers+n_monomers,:]    #substrate concentrations         
# E_tseries =  sol[1+n_polymers+n_monomers:n_polymers+n_monomers+n_microbes,:]'    #reserve biomass mM, when indexing a range seems need to transpose
# V_tseries =  sol[1+n_polymers+n_monomers+n_microbes:n_polymers+n_monomers+2*n_microbes,:]'    #structural biomass mM
# X_tseries =  sol[n_polymers+n_monomers+2*n_microbes+1,:]    #enzyme concentration mM
# CO2_tseries =  sol[n_polymers+n_monomers+2*n_microbes+2,:]  # CO2 cummulative mM
# Bio = E_tseries.+V_tseries
# # N_cells_tseries  = @. Bio[1]*1e-6*12.011/(initb["rhoB"]*initb["Md"])[1] #rhoB is ρ_bulk (1 g/cm^3) and Md is dry_mass, see isolates_batch_init.jl
# N_cells_tseries  = @. Bio.*1e-6.*12.011./(initb["rhoB"].*initb["Md"])[id_isolate]
# BGE_tseries= Bio./(Bio.+CO2_tseries) #Bacterial Growth Efficiecy aka CUE
# # r    = [DEBmicroTrait.growth!(0.0*ones(1), p.metabolism_pars, E_tseries[:][i], V_tseries[:][i]) for i in 1:size(sol.t,1)]
# # r_tseries=0.0*ones(D_tseries)
# # for k in 1:length(sol.t)
# #     r_tseries[k] = r[k]
# # end
# #r    = [DEBmicroTrait.growth!(0.0*ones(1), p.metabolism_pars, [sol[i][2]], [sol[i][3]])[1] for i in 1:size(sol.t,1)]

# #r    = [DEBmicroTrait.growth!(0.0*ones(1), p.metabolism_pars, [sol[i][2]], [sol[i][3]])[1] for i in 1:size(sol.t,1)]

# r = [DEBmicroTrait.growth!(0.0*ones(1), p.metabolism_pars, [E_tseries[i]], [V_tseries[i]])[1] for i in 1:size(sol.t,1)]

# ############## Plotting
# l = @layout [a b c;d e f] #initialize subplot layout

# p1=plot(sol.t, D_tseries',legend=false)#,label=df_metabolites.Formula')
# ylabel!("[Substrate] (mM)")
# xlabel!("Time (hr)")
# xlims!(0,50)

# p2=plot(sol.t,E_tseries, legend=false)
# ylabel!("Reserve (mM)")
# xlabel!("Time (hr)")
# xlims!(0,50)

# p3=plot(sol.t, V_tseries, legend=false)
# ylabel!("[Structural] (mM)")
# xlabel!("Time (hr)")
# xlims!(0,50)

# p4=plot(sol.t, X_tseries, legend=false)
# ylabel!("[Enzyme] (mM)")
# xlabel!("Time (hr)")
# xlims!(0,50)

# p5=plot(sol.t,CO2_tseries, legend=false)
# ylabel!("Cummulative [CO2] (mM)")
# xlabel!("Time (hr)")
# xlims!(0,50)

# p6=plot(sol.t,N_cells_tseries, yscale=:log10,legend=false)
# ylabel!("Number of Cells (mM)")
# xlabel!("Time (hr)")
# xlims!(0,50)
# using Plots.PlotMeasures
# p7=plot(p1, p2, p3, p4, p5, p6, layout = l, size=(1200,800),left_margin=[20mm 0mm])

# l2 = @layout [a b] #initialize subplot layout

# p8=plot(sol.t,BGE_tseries, legend=false)
# ylabel!("BGE")
# xlabel!("Time (hr)")
# xlims!(0,50)

# p9=plot(sol.t, r, legend=false)
# ylabel!("growth rate [1/hr]")
# xlabel!("Time (hr)")
# xlims!(0,50)

# p10= plot(p8, p9, layout = l2, size=(1200,800),left_margin=[20mm 0mm])
# display(p7)
# display(p10)

# # #Added from rhizoshpere_batch_model to calculate growth rate and BGE

# # BGE_tseries       = zeros(n_isolates, n_monomers, n_t) # Bacterial growth efficiency, BGE = BP/ (BP+BR)
# # BR_tseries        = zeros(n_isolates, n_monomers, n_t) # Biomass respiration = assilation respiration + growth respiration + enzyme production respiration + maintenance respiration
# # BP_tseries        = zeros(n_isolates, n_monomers, n_t) # Biomass production = E + V production
# # r_tseries         = zeros(n_isolates, n_monomers, n_t) # realized growth rate [1/h]
# # x_tseries         = zeros(n_isolates, n_monomers, n_t) # constitutive enzyme production rate
# # rG_CO2_tseries    = zeros(n_isolates, n_monomers, n_t) # growth respiration
# # rM_CO2_tseries    = zeros(n_isolates, n_monomers, n_t) # maintenance respiration
# # rX_CO2_tseries    = zeros(n_isolates, n_monomers, n_t) # enzyme production respiration
# # J_EX_tseries      = zeros(n_isolates, n_monomers, n_t) # total enzyme production rate i.e. x*V
# # J_DE_tseries      = zeros(n_isolates, n_monomers, n_t) # assimilation flux
# # J_DE_CO2_tseries  = zeros(n_isolates, n_monomers, n_t) # assimilation respiration
# # J_D_tseries       = zeros(n_isolates, n_monomers, n_t) # uptake rate
# # J_ED_tseries      = zeros(n_isolates, n_monomers, n_t) # recycling of reserve biomass E to substrates D
# # J_V_tseries       = zeros(n_isolates, n_monomers, n_t) # structural biomass turnover rate
# # J_E_tseries       = zeros(n_isolates, n_monomers, n_t) # reserve biomass turnover rate
# # t_tseries         = zeros(n_isolates, n_monomers, n_t) # time (min?)
# # D_tseries         = zeros(n_isolates, n_monomers, n_t) # Substrate concentration
# # E_tseries         = zeros(n_isolates, n_monomers, n_t) # Reserve concentration
# # V_tseries         = zeros(n_isolates, n_monomers, n_t) # Structural biomass concentration
# # X_tseries         = zeros(n_isolates, n_monomers, n_t) # Enzyme concentration
# # CO2_tseries       = zeros(n_isolates, n_monomers, n_t) # CO2 concentration
# # N_cells_tseries   = zeros(n_isolates, n_monomers, n_t) # Conversion from total biomass to number of cells
# # maintenance_tseries    = zeros(n_isolates, n_monomers, n_t) #maintenance respiration

# # for i in 1:n_isolates
# #     for j in 1:n_monomers
# #         for k in 1:length(sol.t)
# #             t_tseries[i,j,k] = sol.t[k]               #time
# #             D_tseries[i,j,k] =  sol[k][1]             
# #             E_tseries[i,j,k] =  sol[k][2]
# #             V_tseries[i,j,k] =  sol[k][3]
# #             X_tseries[i,j,k] =  sol[k][4]
# #             CO2_tseries[i,j,k] =  sol[k][5]
# #             Bio = sol[k][2].+sol[k][3]
# #             N_cells_tseries[i,j,k]  = @. Bio[1]*1e-6*12.011/(initb["rhoB"]*initb["Md"])[1]
# #         end
# #         du   = zeros(p.setup_pars.dim)
# #         BR   = [DEBmicroTrait.batch_model!(du, sol.u[i], p, 0)[end] for i in 1:size(sol.t,1)]
# #         BP   = [DEBmicroTrait.batch_model!(du, sol.u[i], p, 0)[2] + DEBmicroTrait.batch_model!(du, sol.u[i], p, 0)[3] for i in 1:size(sol.t,1)]
# #         BGE  = @. BP/(BP + BR)
# #         for k in 1:length(sol.t)
# #             BGE_tseries[i,j,k] = BGE[k]               #Bacterial growth efficiency aka CUE
# #         end
# #         for k in 1:length(sol.t)
# #             BP_tseries[i,j,k] = BP[k]                 #Biomass production
# #         end
# #         for k in 1:length(sol.t)
# #             BR_tseries[i,j,k] = BR[k]                 #Biomass respiration
# #         end
# #         r    = [DEBmicroTrait.growth!(0.0*ones(1), p.metabolism_pars, [sol[i][2]], [sol[i][3]])[1] for i in 1:size(sol.t,1)]
# #         for k in 1:length(sol.t)
# #             r_tseries[i,j,k] = r[k]
# #         end

# #     end
# # end

# #         for k in 1:length(sol.t)
# #             x, rG_CO2, rM_CO2, rX_CO2 = DEBmicroTrait.growth_production!(r[k], p.metabolism_pars, [sol[k][2]], [sol[k][3]])
# #             x_tseries[i,j,k] = x[1]
# #             rG_CO2_tseries[i,j,k] = rG_CO2[1]
# #             rM_CO2_tseries[i,j,k] = rM_CO2[1]
# #             rX_CO2_tseries[i,j,k] = rX_CO2[1]
# #             J_EX   = DEBmicroTrait.enzyme_production!(x[1], p.metabolism_pars, [sol[k][3]])
# #             J_EX_tseries[i,j,k] = J_EX[1]
# #         end

# #         for k in 1:length(sol.t)
# #             J_DE  = DEBmicroTrait.assimilation!(zeros(1), p.assimilation_pars, [sol[k][1]], [sol[k][3]])
# #             J_DE_tseries[i,j,k] = J_DE[1]
# #             J_DE_CO2 = DEBmicroTrait.assimilation_production!(zeros(1), p.assimilation_pars, [sol[k][1]], [sol[k][3]])
# #             J_DE_CO2_tseries[i,j,k] = J_DE_CO2[1]
# #             J_D      = DEBmicroTrait.uptake!(zeros(1), p.assimilation_pars, [sol[k][1]], [sol[k][3]])
# #             J_D_tseries[i,j,k] = J_D[1]
# #         end

# #         for k in 1:length(sol.t)
# #             J_ED         = DEBmicroTrait.reserve_recycling!(zeros(1), p.turnover_pars, [sol[k][2]])
# #             J_ED_tseries[i,j,k] = J_ED[1]
# #             J_V          = DEBmicroTrait.biomass_turnover!(zeros(1), p.turnover_pars, [sol[k][3]])
# #             J_V_tseries[i,j,k] = J_V[1]
# #             J_E          = DEBmicroTrait.biomass_turnover!(zeros(1), p.turnover_pars, [sol[k][2]])
# #             J_E_tseries[i,j,k] = J_E[1]
# #         end

# #         for k in 1:length(sol.t)
# #             x, rG_CO2, rM_CO2, rX_CO2 = DEBmicroTrait.growth_production!(r[k], p.metabolism_pars, [sol[k][2]], [sol[k][3]])
# #             maintenance_tseries[i,j,k] = rM_CO2[1]./(rG_CO2[1]+rM_CO2[1]+rX_CO2[1])
# #         end
# #     end    
# # end

# # df_out_tseries         = DataFrame()
# # df_out_tseries.time    = vec(t_tseries)
# # df_out_tseries.N_cells = vec(N_cells_tseries)
# # df_out_tseries.BR      = vec(BR_tseries)
# # df_out_tseries.CO2     = vec(CO2_tseries)
# # df_out_tseries.BGE     = vec(BGE_tseries)


# # # #plot N_cells_tseries (number of cells), BR_tseries (respiration), CO2_tseries, and BGE_tseries (Growth efficiency), 
# # #Note need to update file name for different strains

# # CSV.write(joinpath(dir, "files/output2/1_10_R2A_batch_model_tseries_HA54.csv"), df_out_tseries)
