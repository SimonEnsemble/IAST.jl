#=
the Temkin approx
=#
# params
struct TemkinApproxModel<:AdsorptionIsothermModel
	M::Real
	K::Real
    θ::Real
end
TemkinApproxModel(;M::Real=NaN, K::Real=NaN, θ::Real=NaN) = TemkinApproxModel(M, K, θ)

# amount adsorbed as a function of pressure p
function loading(p::Float64, tm::TemkinApproxModel)
    ψ = tm.K * p / (1 + tm.K * p)
    return tm.M * ϕ * (1 + tm.θ * ψ * (ψ - 1))
end

# grand potential (::"spreading pressure") as a function of p
function grand_pot(p, tm::TemkinApproxModel)
    ϕ = 1 + tm.K * p
    return tm.M * (log(ϕ) + tm.θ * (1 + 2 * tm.K * p) / (2 * ϕ^2))
end

function _default_θ_guess(model::TemkinApproxModel,
	                      data::DataFrame,
	                      p_key::String,
	                      l_key::String)
    lm_lang = _default_θ_guess(LangmuirModel(), data, p_key, l_key)
    return TemkinApproxModel(M=lm_lang.M, K=lm_lang.K, θ=0.0)
end
