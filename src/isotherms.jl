abstract type AdsIsoTModel
end

Broadcast.broadcastable(m::AdsIsoTModel) = (m,) # allow loading.(p, model)

function model_to_θ(model::AdsIsoTModel)
	fields = fieldnames(typeof(model))
	return [getfield(model, f) for f in fields]
end

nb_params(model::AdsIsoTModel) = length(model_to_θ(model))

#=
list of isotherm models
=#
include("models/langmuir.jl")
include("models/temkin.jl")
include("models/quadratic.jl")
include("models/dual_site_langmuir.jl")

model_list = [LangmuirModel(), TemkinApproxModel(), QuadraticModel(), DualSiteLangmuirModel()]
