using DEBmicroTrait
using CSV, DataFrames, Statistics
using JLD

########################################
# I/O
dir                     = "DEBSCRIPTS" in keys(ENV) ? ENV["DEBSCRIPTS"] : pwd()
df_isolates             = CSV.read(joinpath(dir, "files/input/isolates2traits.csv"), DataFrame, missingstring="N/A")

########################################
Genome_size     = convert(Array{Float64,1}, df_isolates.Genome_size)
V_cell          = DEBmicroTrait.genome_size_to_cell_volume(Genome_size)
Min_gen_time    = df_isolates.Min_gen_time
gmax            = log(2)./Min_gen_time
Gram_stain      = convert(Array{String,1}, df_isolates.gram_stain)########################################

########################################
dry_mass        = 0.47*DEBmicroTrait.cell_volume_to_dry_mass(V_cell, gmax, Gram_stain)
ρ_bulk          = 1.0 # g/cm^3
N_cells         = 1e2 #Where N cells is number of cells per gram of soil/ solution
Bio_0           = N_cells*1e6*ρ_bulk*dry_mass./12.011
########################################

########################################
# I/O, change to my directory
save("./files/output/isolates_batch_init.jld", "Bio0", Bio_0, "Md", dry_mass, "rhoB", ρ_bulk, "Ncells", N_cells)
########################################
