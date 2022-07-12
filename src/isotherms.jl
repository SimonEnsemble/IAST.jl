abstract type AdsorptionIsothermModel
end

#=
the Langmuir adsorption isotherm model
=#
# params
struct LangmuirModel<:AdsorptionIsothermModel
	M::Float64
	K::Float64
end

# amount adsorbed
n(p::Float64, lm::LangmuirModel) = lm.M * lm.K * p / (1 + lm.K * p)

# spreading pressure
π(p::Float64, lm::LangmuirModel) = lm.M * log(1.0 + lm.K * p)
