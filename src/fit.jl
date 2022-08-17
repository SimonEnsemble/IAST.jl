function identify_params(ads_data::AdsIsoTData, model::AdsIsoTModel)
	θ₀ = model_to_θ(model)
	if any(isnan.(θ₀))
        model_guess = _default_θ_guess(ads_data, model)
		θ₀ = model_to_θ(model_guess)
	end
	
	# construct sum of square errors loss function
	function ℓ(θ)
		# construct model
		model = typeof(model)(θ...)
		
		loss = 0.0
		for row in eachrow(ads_data.data)
			# unpack data point
			pᵢ, nᵢ = row[ads_data.p_key], row[ads_data.l_key]
			
			# predicted loading
			n̂ᵢ = loading(pᵢ, model)

			# increment loss.
			loss += (nᵢ - n̂ᵢ)^2
		end
		return loss
	end

	# minimize loss
	res = optimize(ℓ, θ₀)

	# return model
	return typeof(model)(res.minimizer...)
end
