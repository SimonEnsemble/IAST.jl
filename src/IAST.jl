module IAST

using NLsolve, Statistics, LinearAlgebra, DataFrames, Optim
include("data.jl")
include("isotherms.jl")
include("fit.jl")
include("calcs.jl")

export AdsIsoTData, pressures, # data.jl
       # isotherms.jl
       loading, grand_pot, AdsIsoTModel,
       LangmuirModel, TemkinApproxModel, QuadraticModel, DualSiteLangmuirModel,
       # calcs.jl
       iast, 
       # fit.jl
       identify_params
