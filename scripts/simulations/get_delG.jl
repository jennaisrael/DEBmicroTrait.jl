using DEBmicroTrait
using CSV, DataFrames, Statistics
using JLD
using DifferentialEquations

dir                     = "DEBSCRIPTS" in keys(ENV) ? ENV["DEBSCRIPTS"] : pwd()
df_isolates             = CSV.read(joinpath(dir, "files/input/isolates2traits.csv"), DataFrame, missingstring="N/A")
df_metabolites          = CSV.read(joinpath(dir, "files/input/root_exudates.csv"), DataFrame, missingstring="N/A")



# thermodynamic efficiency
dGcox              = zeros(83)
dGcat              = zeros(83)
dGAn               = zeros(83)
λ_base             = zeros(83)
N_C                = zeros(83)
eta                = zeros(83)
chemFormBiom       = [1, 1.64, 0.2, 0.38, 0, 0, 0]

for i in 1:83
    elementstring = convert(String,df_metabolites.Formula[i])
    N_C[i] = DEBmicroTrait.extract_composition(elementstring)[1]
    out = DEBmicroTrait.get_lambda(elementstring, chemFormBiom)
    dGcox[i] = out[2][3]
    dGcat[i] = out[2][5]
    dGAn[i]  = out[2][8]
    λ_base[i]     = out[1][1]  # where lambda is how many times the catabolic reaction 
    #needs to run to provide the energy required for the synthsisi of a unit C-mole of biomass (Song et al. 2020)
    eta      = @. dGAn/(λ_base*dGcat)
end

df_all = DataFrame()

monomer = Array{String}(undef,39,83)
for i in 1:83
     monomer[:,i] .= df_metabolites.Name[i]
 end
df_all.monomer = vec(monomer)
ontology = Array{String}(undef,39,83)
for i in 1:83
     ontology[:,i] .= df_metabolites.Ontology[i]
 end
df_all.ontology = vec(ontology)
eta_eff = zeros(39,83)
for i in 1:83
     eta_eff[:,i] .= eta[i]
 end
df_all.eta = vec(eta_eff)

lambda_eff = zeros(39,83)
for i in 1:83
     lambda_eff[:,i] .=  λ_base[i]
 end
df_all.lambda = vec(lambda_eff)



CSV.write(joinpath(dir, "files/output2/substrates_deltaG.csv"), df_all)
