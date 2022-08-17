#=
the Temkin approx
=#
# params
struct TemkinApproxModel<:AdsIsoTModel
    M::Real
    K::Real
    θ::Real
end
TemkinApproxModel(;M::Real=NaN, K::Real=NaN, θ::Real=NaN) = TemkinApproxModel(M, K, θ)

# amount adsorbed as a function of pressure p
function loading(p::Real, tm::TemkinApproxModel)
    ψ = tm.K * p / (1 + tm.K * p)
    return tm.M * ψ * (1 + tm.θ * ψ * (ψ - 1))
end

# grand potential (::"spreading pressure") as a function of p
function grand_pot(p::Real, tm::TemkinApproxModel)
    ϕ = 1 + tm.K * p
    return tm.M * (log(ϕ) + tm.θ * (1 + 2 * tm.K * p) / (2 * ϕ^2))
end

function _default_θ_guess(ads_data::AdsIsoTData, model::TemkinApproxModel)
    lm_lang = _default_θ_guess(ads_data, LangmuirModel())
    return TemkinApproxModel(M=lm_lang.M, K=lm_lang.K, θ=0.0)
end
