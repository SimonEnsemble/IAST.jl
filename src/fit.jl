function fit(data::DataFrame, 
	         p_key::String, 
	         l_key::String,
	         model::AdsorptionIsothermModel)
	θ₀ = model_to_θ(model)
	if any(isnan.(θ₀))
		model_guess = _default_θ_guess(model, data, p_key, l_key)
		θ₀ = model_to_θ(model_guess)
	end
	
	# construct sum of square errors loss function
	function ℓ(θ)
		# construct model
		model = typeof(model)(θ...)
		
		loss = 0.0
		for row in eachrow(data)
			# unpack data point
			pᵢ, nᵢ = row[p_key], row[l_key]
			
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
