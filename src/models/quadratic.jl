#=
the quadratic isotherm model
=#
# params
struct QuadraticModel<:AdsIsoTModel
    M::Real
    K::Real
    ϕ::Real
end
QuadraticModel(;M::Real=NaN, K::Real=NaN, ϕ::Real=NaN) = QuadraticModel(M, K, ϕ)

# amount adsorbed as a function of pressure p
function loading(p::Real, qm::QuadraticModel)
    Kp = qm.K * p
    return 2 * qm.M * Kp * (1 + qm.ϕ * Kp) / (1 + 2 * Kp + qm.ϕ * Kp^2) 
end

# grand potential (::"spreading pressure") as a function of p
function grand_pot(p::Real, qm::QuadraticModel)
    return qm.M * log(1 + 2 * qm.K * p + qm.ϕ * (qm.K * p)^2)
end

function _default_θ_guess(ads_data::AdsIsoTData, model::QuadraticModel)
    lm = _default_θ_guess(ads_data, LangmuirModel())
    return QuadraticModel(M=lm.M / 2, K=lm.K, ϕ=1.0)
end
