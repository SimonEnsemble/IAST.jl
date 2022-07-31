module IAST

using NLsolve, Statistics, LinearAlgebra, DataFrames, Optim
include("isotherms.jl")
include("fit.jl")
include("calcs.jl")

export loading, grand_pot, AdsorptionIsothermModel, LangmuirModel, TemkinApproxModel, iast, identify_params

end
