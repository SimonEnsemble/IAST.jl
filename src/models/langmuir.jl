#=
the Langmuir adsorption isotherm model
=#
# params
struct LangmuirModel<:AdsIsoTModel
    M::Real
    K::Real
end
LangmuirModel(;M::Real=NaN, K::Real=NaN) = LangmuirModel(M, K)

# amount adsorbed as a function of pressure p
loading(p::Real, lm::LangmuirModel) = lm.M * lm.K * p / (1 + lm.K * p)

# grand potential (::"spreading pressure") as a function of p
grand_pot(p::Real, lm::LangmuirModel) = lm.M * log(1.0 + lm.K * p)

function _default_θ_guess(ads_data::AdsIsoTData, model::LangmuirModel)
    # use first two data points to get the slope
    H = mean([ads_data.data[i, ads_data.l_key] / ads_data.data[i, ads_data.p_key] for i = 1:2])

    # estimate M as 10% more than max observed loading
    M = maximum(ads_data.data[:, ads_data.l_key])
    K = H / M

    return LangmuirModel(K=K, M=M)
end
