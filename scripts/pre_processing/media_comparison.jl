using DEBmicroTrait
using CSV, DataFrames, Statistics
using Roots
using JLD

########################################
# I/O
dir                     = "DEBSCRIPTS" in keys(ENV) ? ENV["DEBSCRIPTS"] : pwd()

#read in the media files with only unique molecules
df_R2A_1_10              = CSV.read(joinpath(dir, "files/output2/1_100_R2A_medium_unique.csv"), DataFrame, missingstring="N/A")
df_LB                    = CSV.read(joinpath(dir, "files/output2/LB_medium_unique.csv"), DataFrame, missingstring="N/A")

#left out the minerals, add the P containing ones
# LB is Casein Peptone 10g/L, Yeast extract 5g/L, Sodium Chloride 10g/L, (https://www.fishersci.com/shop/products/lb-broth-miller-granulated/p-4245497)
push!(df_LB, ("PO4", "phosphate", "P Source", 94.971, 4.659)) # 2.6956 mmol/L from casein hydrolysate and 1.964 from yeast extract

#R2A has P in Casein Peptone, yeast extract, and Proteose Peptone
push!(df_R2A_1_10, ("PO4", "phosphate", "P Source", 94.971, 0.03648)) #see calculations in 1_10_R2A_media_with_P.csv

#get the stoichiometry of the media,
df_LB.Formula  = convert.(String, df_LB.Formula)



N_CNP_LB                     = zeros(size(df_LB.Compound,1),3)
for i in 1:size(df_LB.Compound,1)
    elementstring       = df_LB.Formula[i]
    #Carbon is index 1, Nitrogen is 3, and P is 6 (see extract_composition function in thermostoichwizard.jl)
    N_CNP_LB[i,1]              = DEBmicroTrait.extract_composition(elementstring)[1]
    N_CNP_LB[i,2]              = DEBmicroTrait.extract_composition(elementstring)[3] 
    N_CNP_LB[i,3]              = DEBmicroTrait.extract_composition(elementstring)[6]
end
#Scale each row by the concentration of the metabolite
C_CNP_LB=N_CNP_LB.*df_LB.Concentration_sum
#Sum the rows
LB_ratio_raw=sum(C_CNP_LB; dims=1)
LB_ratio=LB_ratio_raw./minimum(LB_ratio_raw)


#get the stoichiometry of the media


df_R2A_1_10.Formula  = convert.(String, df_R2A_1_10.Formula)
N_CNP_R2A                  = zeros(size(df_R2A_1_10.Compound,1),3)
for i in 1:size(df_R2A_1_10.Compound,1)
    elementstring       = df_R2A_1_10.Formula[i]
    #Carbon is index 1, Nitrogen is 3, and P is 6 (see extract_composition function in thermostoichwizard.jl)
    N_CNP_R2A[i,1]              = DEBmicroTrait.extract_composition(elementstring)[1]
    N_CNP_R2A[i,2]              = DEBmicroTrait.extract_composition(elementstring)[3] 
    N_CNP_R2A[i,3]              = DEBmicroTrait.extract_composition(elementstring)[6]
end

#Scale each row by the concentration of the metabolite
C_CNP_R2A=N_CNP_R2A.*df_R2A_1_10.Concentration_sum
#Sum the rows
R2A_ratio_raw=sum(C_CNP_R2A; dims=1)
R2A_ratio=R2A_ratio_raw./minimum(R2A_ratio_raw)
