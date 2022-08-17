#=
the dual-site Langmuir adsorption isotherm model
=#
# params
struct DualSiteLangmuirModel<:AdsorptionIsothermModel
    M₁::Real
    K₁::Real
    M₂::Real
    K₂::Real
end
DualSiteLangmuirModel(;M₁::Real=NaN, K₁::Real=NaN, M₂::Real=NaN, K₂::Real=NaN) = DualSiteLangmuirModel(M₁, K₁, M₂, K₂)

# amount adsorbed as a function of pressure p
loading(p::Real, m::DualSiteLangmuirModel) = m.M₁ * m.K₁ * p / (1 + m.K₁ * p) + m.M₂ * m.K₂ * p / (1 + m.K₂ * p)

# grand potential (::"spreading pressure") as a function of p
grand_pot(p::Real, m::DualSiteLangmuirModel) = m.M₁ * log(1.0 + m.K₁ * p) + m.M₂ * log(1.0 + m.K₂ * p)

function _default_θ_guess(ads_data::AdsorptionIsothermData, model::DualSiteLangmuirModel)
    lm = _default_θ_guess(ads_data, LangmuirModel())

    return DualSiteLangmuirModel(K₁=lm.K * 2, M₁=lm.M / 2, K₂=lm.K / 2, M₂=lm.M / 2)
end
