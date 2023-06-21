#Create a master file to change V_cell and generate necessary output files
#NOTE: currently configured to run one cell at a time

using DEBmicroTrait
using CSV, DataFrames, Statistics
using Roots
using JLD

#V_cell                  = DEBmicroTrait.genome_size_to_cell_volume(convert(Array{Float64,1}, df_isolates.Genome_size))
V_cell                  = 2.99* 10^(-18)*ones(1)#this needs to be a vector, m^3
fileend                 = "LB_HB15_live_mean.jld"

#####Assimilation
########################################
# I/O
dir                     = "DEBSCRIPTS" in keys(ENV) ? ENV["DEBSCRIPTS"] : pwd()
df_isolates_all            = CSV.read(joinpath(dir, "files/input/isolates2traits.csv"), DataFrame, missingstring="N/A")
###########Add a step here that aggregates the duplicate entries in the media input file so the "Formula" column entries are unique
###########and the sum of all the different sources
df_media_raw            = CSV.read(joinpath(dir, "files/input/LB_medium.csv"), DataFrame, missingstring="N/A")
gdf = groupby(df_media_raw, [:Formula, :Compound, :Ontology, :Molecular_weight]) #need columns named "Formula", "Name"(Name should probably be "Compouund" instead), "Ontology", "Molecular_weight" and "Concentration" (used in the batch model)
unique_metabolites = combine(gdf, :Concentration => sum)
#write this to a CSV
CSV.write(joinpath(dir, "files/output2/LB_medium_unique.csv"),unique_metabolites)
df_metabolites_all          = CSV.read(joinpath(dir, "files/output2/LB_medium_unique.csv"), DataFrame, missingstring="N/A")
########################################
# metabolite traits, replace "Name" with "Compound"
df_metabolites_all.Formula  = convert.(String, df_metabolites_all.Formula)
N_C                     = zeros(size(df_metabolites_all.Compound,1))
for i in 1:size(df_metabolites_all.Compound,1)
    elementstring       = df_metabolites_all.Formula[i]
    N_C[i]              = DEBmicroTrait.extract_composition(elementstring)[1]
end
########################################
###EDIT Trim df_isolates_all and df_metabolites_all to be just one strain and one substrate (HA54 on glycine)
n_isolates=1
#df_metabolites=filter(:Compound => ==("serine"), df_metabolites_all)
df_metabolites=df_metabolites_all # this time run for all LB 

df_isolates=filter(:Abbreviation => ==("HB15"), df_isolates_all)

# isolate traits
# V_cell                  = DEBmicroTrait.genome_size_to_cell_volume(convert(Array{Float64,1}, df_isolates.Genome_size))
Min_gen_time            = df_isolates.Min_gen_time
Gram_stain              = convert(Array{String,1}, df_isolates.gram_stain)
rrn_copies              = convert(Array{Float64,1}, df_isolates.rRNA_genes)
y_EM                    = ones(size(V_cell,1))

z_sugars                = reshape(convert(Array{Float64,1}, df_isolates.z_sugars./df_isolates.Genome_size*1e6),1,n_isolates)
z_organics              = reshape(convert(Array{Float64,1}, df_isolates.z_organic_acids./df_isolates.Genome_size*1e6),1,n_isolates)
z_aminos                = reshape(convert(Array{Float64,1}, df_isolates.z_amino_acids./df_isolates.Genome_size*1e6), 1,n_isolates)
z_fattys                = reshape(convert(Array{Float64,1}, df_isolates.z_fatty_acids./df_isolates.Genome_size*1e6),1,n_isolates)
z_nucleos               = reshape(convert(Array{Float64,1}, df_isolates.z_nucleotides./df_isolates.Genome_size*1e6),1,n_isolates)
z_auxins                = reshape(convert(Array{Float64,1}, df_isolates.z_auxins./df_isolates.Genome_size*1e6),1,n_isolates)
genome_distr            = vcat(z_sugars, z_organics, z_aminos, z_fattys, z_nucleos, z_auxins)
########################################
# estimate transporter density
ρ_ps                    = zeros(size(df_metabolites.Compound,1), size(V_cell,1))
y_DEs                   = zeros(size(df_metabolites.Compound,1), size(V_cell,1))



for j in 1:size(df_metabolites.Compound,1)
    if df_metabolites.Ontology[j] == "Sugars"
        for i in 1:size(V_cell,1)
            find_ρ(x)   = DEBmicroTrait.constrain_transporter_density_cost(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], df_metabolites.Formula[j])
            ρ_p         = Roots.find_zero(find_ρ, 1.0)
            closure     = genome_distr[:,i]./sum(genome_distr[:,i])
            ρ_ps[j,i]   = ρ_p[1].*closure[1]
            y_DE        = DEBmicroTrait.yield_transporter_density_cost(ρ_ps[j,i], [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], df_metabolites.Formula[j])
            y_DEs[j,i]  = y_DE[1]
        end
    elseif df_metabolites.Ontology[j] == "Organic acids"
        for i in 1:size(V_cell,1)
            find_ρ(x)   = DEBmicroTrait.constrain_transporter_density_cost(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], df_metabolites.Formula[j])
            ρ_p         = Roots.find_zero(find_ρ, 1.0)
            closure     = genome_distr[:,i]./sum(genome_distr[:,i])
            ρ_ps[j,i]   = ρ_p[1].*closure[2]
            y_DE        = DEBmicroTrait.yield_transporter_density_cost(ρ_ps[j,i], [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], df_metabolites.Formula[j])
            y_DEs[j,i]  = y_DE[1]
        end
    elseif df_metabolites.Ontology[j] == "Amino acids"
        for i in 1:size(V_cell,1)
            find_ρ(x)   = DEBmicroTrait.constrain_transporter_density_cost(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], df_metabolites.Formula[j])
            ρ_p         = Roots.find_zero(find_ρ, 1.0)
            closure     = genome_distr[:,i]./sum(genome_distr[:,i])
            ρ_ps[j,i]   = ρ_p[1].*closure[3]
            y_DE        = DEBmicroTrait.yield_transporter_density_cost(ρ_ps[j,i], [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], df_metabolites.Formula[j])
            y_DEs[j,i]  = y_DE[1]
        end
    elseif df_metabolites.Ontology[j] == "Fatty acids"
        for i in 1:size(V_cell,1)
            find_ρ(x)   = DEBmicroTrait.constrain_transporter_density_cost(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], df_metabolites.Formula[j])
            ρ_p         = Roots.find_zero(find_ρ, 1.0)
            closure     = genome_distr[:,i]./sum(genome_distr[:,i])
            ρ_ps[j,i]   = ρ_p[1].*closure[4]
            y_DE        = DEBmicroTrait.yield_transporter_density_cost(ρ_ps[j,i], [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], df_metabolites.Formula[j])
            y_DEs[j,i]  = y_DE[1]
        end
    elseif df_metabolites.Ontology[j] == "Nucleotides"
        for i in 1:size(V_cell,1)
            find_ρ(x)   = DEBmicroTrait.constrain_transporter_density_cost(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], df_metabolites.Formula[j])
            ρ_p         = Roots.find_zero(find_ρ, 1.0)
            closure     = genome_distr[:,i]./sum(genome_distr[:,i])
            ρ_ps[j,i]   = ρ_p[1].*closure[5]
            y_DE        = DEBmicroTrait.yield_transporter_density_cost(ρ_ps[j,i], [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], df_metabolites.Formula[j])
            y_DEs[j,i]  = y_DE[1]
        end
    else
        for i in 1:size(V_cell,1)
            find_ρ(x)   = DEBmicroTrait.constrain_transporter_density_cost(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], df_metabolites.Formula[j])
            ρ_p         = Roots.find_zero(find_ρ, 1.0)
            closure     = genome_distr[:,i]./sum(genome_distr[:,i])
            ρ_ps[j,i]   = ρ_p[1].*closure[6]
            y_DE        = DEBmicroTrait.yield_transporter_density_cost(ρ_ps[j,i], [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], df_metabolites.Formula[j])
            y_DEs[j,i]  = y_DE[1]
        end
    end
end

ρ_ps[ρ_ps.==0.0] .= 1e-12
median(ρ_ps)

########################################
# calculate uptake traits
N_SB              = DEBmicroTrait.transporter_density_to_monomer_uptake_sites(V_cell, ρ_ps, Min_gen_time, Gram_stain)
Vmax              = @. 180.0*60^2*N_SB.*N_C
#
D_S               = DEBmicroTrait.aqueous_diffusivity(df_metabolites.Molecular_weight)
K_D               = DEBmicroTrait.specific_reference_affinity(V_cell, ρ_ps, D_S)
#
a_s               = Vmax./K_D

########################################
# I/O
# save("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/output/isolates_assimilation_10_R2A.jld", "rho", ρ_ps, "NSB", N_SB, "KD", K_D, "yEM", y_EM, "yDE", y_DEs, "NC", N_C)
save(string("./files/output/isolates_assimilation_",fileend), "rho", ρ_ps, "NSB", N_SB, "KD", K_D, "yEM", y_EM, "yDE", y_DEs, "NC", N_C)


####### batch_init########################################################################################################################
########################################
#Genome_size     = convert(Array{Float64,1}, df_isolates.Genome_size)
#V_cell          = DEBmicroTrait.genome_size_to_cell_volume(Genome_size)
Min_gen_time    = df_isolates.Min_gen_time
gmax            = log(2)./Min_gen_time
Gram_stain      = convert(Array{String,1}, df_isolates.gram_stain)

########################################
dry_mass        = 0.47*DEBmicroTrait.cell_volume_to_dry_mass(V_cell, gmax, Gram_stain)
ρ_bulk          = 1.0 # g/cm^3
N_cells         = 1e2 #Where N cells is number of cells per gram of soil/ solution
Bio_0           = N_cells*1e6*ρ_bulk*dry_mass./12.011
########################################

########################################
# I/O, change to my directory
save(string("./files/output/isolates_batch_init_",fileend), "Bio0", Bio_0, "Md", dry_mass, "rhoB", ρ_bulk, "Ncells", N_cells)
########################################


#######enzymes########################################################################################################################
########################################
#V_cell          = DEBmicroTrait.genome_size_to_cell_volume(convert(Array{Float64,1}, df_isolates.Genome_size))
Min_gen_time    = df_isolates.Min_gen_time
Gram_stain      = convert(Array{String,1}, df_isolates.gram_stain)
zh              =  df_isolates.z_hydrolases./df_isolates.Genome_size*1e6
α_X             =  1e-2*(df_isolates.z_hydrolases./df_isolates.Genome_size*1e6)./maximum(df_isolates.z_hydrolases./df_isolates.Genome_size*1e6)
########################################

########################################
# I/O
save(string("./files/output/isolates_enzymes_",fileend), "zh", zh, "alpha", α_X)

########maintenance########################################################################################################################
########################################
#Genome_size     = convert(Array{Float64,1}, df_isolates.Genome_size)
#V_cell          = DEBmicroTrait.genome_size_to_cell_volume(Genome_size)
rrn_copies      = convert(Array{Float64,1}, df_isolates.rRNA_genes)
Min_gen_time    = df_isolates.Min_gen_time
Gram_stain      = convert(Array{String,1}, df_isolates.gram_stain)
########################################

########################################
k_M    = DEBmicroTrait.cell_volume_to_specific_maintenance_rate(V_cell, Min_gen_time, Gram_stain)
k_M_med = median(k_M)
########################################

########################################
# I/O
save(string("./files/output/isolates_maintenance_",fileend), "kM", k_M)
########protein synthesis########################################################################################################################
########################################
#Genome_size     = convert(Array{Float64,1}, df_isolates.Genome_size)
#V_cell          = DEBmicroTrait.genome_size_to_cell_volume(Genome_size)
rrn_copies      = convert(Array{Float64,1}, df_isolates.rRNA_genes)
Min_gen_time    = df_isolates.Min_gen_time
Gram_stain      = convert(Array{String,1}, df_isolates.gram_stain)
########################################

########################################
gmax            = log(2)./Min_gen_time
V_p             = DEBmicroTrait.cell_volume_to_protein_volume(V_cell)
V_r             = DEBmicroTrait.cell_volume_to_ribosome_volume(V_cell, gmax)
k_E             = DEBmicroTrait.translation_power(V_p, V_r, Min_gen_time)
########################################

########################################
y_EV            = DEBmicroTrait.relative_translation_efficiency_regression(rrn_copies)
########################################

########################################
# I/O
save(string("./files/output/isolates_protein_synthesis_",fileend), "kE", k_E, "yEV", y_EV, "mingt", Min_gen_time)
########################################


######### turnover########################################################################################################################

########################################
#Genome_size     = convert(Array{Float64,1}, df_isolates.Genome_size)
#V_cell          = DEBmicroTrait.genome_size_to_cell_volume(Genome_size)
Min_gen_time    = df_isolates.Min_gen_time
gmax            = log(2)./Min_gen_time
Gram_stain      = convert(Array{String,1}, df_isolates.gram_stain)
########################################

########################################
γ_V_0 = DEBmicroTrait.max_specific_death_rate(Min_gen_time::Vector{Float64})
########################################

########################################
dry_mass        = 0.47*DEBmicroTrait.cell_volume_to_dry_mass(V_cell, gmax, Gram_stain)
ρ_bulk          = 1.21 # g/cm^3
Bio_0           = 1e9*1e6*ρ_bulk*dry_mass./12.011
γ_V_1           = median(Bio_0)
########################################

########################################
# I/O
save(string("./files/output/isolates_turnover_",fileend), "gV0", γ_V_0, "gV1", γ_V_1)
########################################