#Isolate the NOSC calculation part of the thermostoichwizard
using DEBmicroTrait
using CSV, DataFrames, Statistics
using Roots
using JLD
using Plots

function find_nosc(elementstring::String)
    chemical_indices = DEBmicroTrait.extract_composition(elementstring)
    a = chemical_indices[1] #C
    b = chemical_indices[2] #H
    c = chemical_indices[3] #N
    d = chemical_indices[4] #O
    e = chemical_indices[5] #S
    f = chemical_indices[6] #P
    z = 0

  ne = -z+4*a+b-3*c-2*d+5*e-2*f  # number of electrons transferred in D
  nosc = -(ne/a) + 4
  return nosc
  #print(Base.strcat("NOSC=",string(nosc)))
end

dir                     = "DEBSCRIPTS" in keys(ENV) ? ENV["DEBSCRIPTS"] : pwd()
df_metabolites          = CSV.read(joinpath(dir, "files/input/root_exudates.csv"), DataFrame, missingstring="N/A")
df_metabolites.Formula  = convert.(String, df_metabolites.Formula)
df_metabolites.Name     = convert.(String, df_metabolites.Name)
df_metabolites.NOSC     = df_metabolites.Molecular_weight.*0


for j in 1:length(df_metabolites.Formula)
    comp=df_metabolites.Formula[j]
    df_metabolites.NOSC[j] = find_nosc(comp)
end

#filter out the ones with N
df_metabolites.N_flag    = occursin.("N",df_metabolites.Formula)
df_metabolites.N_flag    = convert.(Float64,df_metabolites.N_flag)
#df_metabolites_filtered  = df_metabolites[df_metabolites.N_flag .< 1.0, :]
df_metabolites_filtered = subset(df_metabolites, :N_flag => ByRow(<(1.0)))
df_metabolites_filtered = subset(df_metabolites_filtered, :NOSC => ByRow(!=(0)))

plot(bar(df_metabolites_filtered.NOSC,legend=false, xrotation=45))
plot!(xticks = ([0:1:length(df_metabolites_filtered.Formula);], df_metabolites_filtered.Name))

#el = "C4H9NO2"
#el = "C6H12O6"
#el = "CO2"
#el = "C4H7NO4" # aspartic acid
#el = "CH3OH" #methanol, but gives the wrong answer need to combine same elements into a single number
#el = "CH4O" #methanol rewritten
#el = "C8H8O3" #vanilin
#el = "C10H9NO2" #indole-3 acetic acid
#el = "C7H13NO2" # stachydrine
#el = "C7H13NO3" #bentonicine
#el = "C15H20O4" #absisic acid
#el = "C5H11NO2" #betaine
#el = "C4H9NO3" #homeserine
#el = "C11H12N2O2" #L-tryptophan
#el = "C6H9N3O2" # L-histidine

#chemFormBiom = [1, 1.8, 0.2, 0.5, 0, 0, 0]

 # CHNOSP
#@test chemical_composition == [4,9,1,2,0,0]

#read in root exudate file and add a column that has the NOSC and a 1 if its only CHO

#find_nosc(el)
