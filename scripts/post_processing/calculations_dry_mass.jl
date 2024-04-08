using DEBmicroTrait
using CSV, DataFrames, Statistics
using GLM
using Roots

###########################################################

# # Thiobacillus
# Genome_size = [2.85e6]
# V_cell = DEBmicroTrait.genome_size_to_cell_volume(Genome_size)
# Gram_stain = ["-"]
# rrn_copies = [2.0]
# gmax = DEBmicroTrait.gmax_regression(rrn_copies)
# dw = DEBmicroTrait.cell_volume_to_dry_mass(V_cell, gmax, Gram_stain)*12.01

########################################
# I/O
dir                     = "DEBSCRIPTS" in keys(ENV) ? ENV["DEBSCRIPTS"] : pwd()
df_isolates             = CSV.read(joinpath(dir, "files/input/isolates2traits.csv"), DataFrame, missingstring="N/A")
id_isolate = 7 #HA54 is 7, HB15 is 30
n_isolates = length(id_isolate)
Genome_size= [convert(Float64,df_isolates[id_isolate,"Genome_size"])]
V_cell= DEBmicroTrait.genome_size_to_cell_volume(Genome_size)
Gram_stain=[convert(String,df_isolates[id_isolate,"gram_stain"])]
rrn_copies=[convert(Float64,df_isolates[id_isolate,"rRNA_genes"])]
gmax = DEBmicroTrait.gmax_regression(rrn_copies)
dw = DEBmicroTrait.cell_volume_to_dry_mass(V_cell, gmax, Gram_stain) # in grams
