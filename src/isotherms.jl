abstract type AdsorptionIsothermModel
end

Broadcast.broadcastable(m::AdsorptionIsothermModel) = (m,) # allow loading.(p, model)

function model_to_θ(model::AdsorptionIsothermModel)
	fields = fieldnames(typeof(model))
	return [getfield(model, f) for f in fields]
end

#=
list of isotherm models
=#
include("models/langmuir.jl")
include("models/temkin.jl")
include("models/quadratic.jl")
include("models/dual_site_langmuir.jl")
