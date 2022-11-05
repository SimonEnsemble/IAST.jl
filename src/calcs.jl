# express difference in spreading pressures
# one entry for each component of the gas mixture. inputs are:
# - Δπ: differences in spreading pressures
# - x: mol fractions in adsorbed phase
# - p: partial pressures
# - aim: adsorption isotherm models
function Δπ!(Δπ, x, p, aim)
             #x::Vector{Float64}, 
             #p::Vector{Float64},
             #aim::Vector{<:AdsIsoTModel})
    n_c = length(x)
    @assert length(Δπ) == n_c - 1
    for c = 1:n_c-1
        Δπ[c] = grand_pot(p[c] / x[c], aim[c]) - grand_pot(p[c+1] / x[c+1], aim[c+1])
    end
end

# conduct IAST calculation
# inputs:
# - p: partial pressures in gas phase
# - models: adsorption models
# output:
# - a: vector of loadings in adsorbed phase
function iast(p::Vector{Float64}, models::Vector{<:AdsIsoTModel};
              p_warn::Array{Float64}=[Inf for _ = 1:length(p)])
    n_c = length(p) # number of components
    
    # guess adsorbed phase mol fraction
    # treat as competitive Langmuir
    M̃ = mean([loading(1000.0, models[c]) for c = 1:n_c]) # high pressure loading
    K̃ = [loading(1e-4, models[c]) / 1e-4 for c = 1:n_c] / M̃  # Henry / M = k
    x_guess = [K̃[c] * p[c] / (1 + dot(K̃, p)) for c = 1:n_c]
    x_guess /= sum(x_guess)

    _Δπ!(Δπ::Vector{Float64}, x::Vector{Float64}) = Δπ!(Δπ, vcat(x, [1-sum(x)]), p, models)

    distance_to_0_or_1 = minimum([minimum(1.0 .- x_guess), minimum(x_guess)])
    res = nlsolve(_Δπ!, x_guess[1:end-1], factor=distance_to_0_or_1)#, autodiff=:forward)
    @assert res.f_converged

    x = vcat(res.zero, [1-sum(res.zero)])
    @assert all(x .> 0.0) && all(x .< 1.0)

    p₀ = p ./ x
    for c = 1:n_c
        if p₀[c] > p_warn[c]
            @warn "extrapolating component $c to $(p₀[c])"
        end
    end

    nₜ = 1.0 / sum(x[c] / loading(p₀[c], models[c]) for c = 1:n_c)

    return x * nₜ
end

function iast(partial_pressures::Dict{String, Float64},
              models::Dict{String, <:AdsIsoTModel};
              p_warn::Dict{String, Float64}=Dict{String, Float64}())
    # unpack list of gases to use for ordering
    gases = [gas for gas in keys(partial_pressures)]

    # sanity checks
    @assert length(models) == length(partial_pressures)
    for gas in gases
        @assert haskey(models, gas)
    end
    for gas in keys(p_warn)
        @assert gas in gases
    end

    # handle p warn defaults
    for gas in gases
        if ! haskey(p_warn, gas)
            p_warn[gas] = Inf
        end
    end

    a = iast([partial_pressures[gas] for gas in gases],
             [models[gas] for gas in gases],
             p_warn=[p_warn[gas] for gas in gases])
    return Dict(gas => aᵢ for (gas, aᵢ) in zip(gases, a))
end
