### A Pluto.jl notebook ###
# v0.19.11

using Markdown
using InteractiveUtils

# ╔═╡ 8c4dd1f8-0fa8-11ed-2080-0fa9ccb116fa
begin
	using Pkg
	Pkg.activate()
	push!(LOAD_PATH, "../../IAST.jl/src")
	using CairoMakie, Distributions, PlutoUI, CSV, DataFrames, IAST, ColorSchemes, Printf, StatsBase
end

# ╔═╡ 4e26b473-510c-4fbc-9bc6-283cf88e5f2a
gases = ["methane", "ethane"]

# ╔═╡ b1a58980-60da-4e75-b9b3-46ce0c31cab5
data = Dict(gas => CSV.read("IRMOF-1_$(gas)_isotherm_298K.csv", DataFrame) 
             for gas in gases)

# ╔═╡ ce6bcada-a06f-476b-8255-e651e60dfcec
p_key = "Pressure(bar)"

# ╔═╡ c43ae81a-8686-4237-b08b-06f3020cdd0b
l_key = "Loading(mmol/g)"

# ╔═╡ 8b827180-af14-4590-8f23-13a219f84a72
p_warn = Dict(gas => maximum(data[gas][:, l_key]) for gas in gases)

# ╔═╡ 3b72d972-b4da-480f-beff-2dc1e670f04e
models = Dict(
	"methane" => identify_params(data["methane"], p_key, l_key, LangmuirModel()),
	"ethane"  => identify_params(data["ethane"], p_key, l_key, LangmuirModel())
)

# ╔═╡ ad050898-0d2f-4a6f-9c75-c4d093b16e04
identify_params(data["methane"], p_key, l_key, LangmuirModel())

# ╔═╡ 5b1cb128-4714-4024-bcc7-5fb4af9dff52
begin
	fig = Figure()
	ax  = Axis(fig[1, 1], xlabel=p_key, ylabel=l_key)
	for gas in gases
		scatter!(data[gas][:, p_key], data[gas][:, l_key], label=gas)

		ps = range(0, 150, length=100)
		lines!(ps, loading.(ps, models[gas]))
	end
	axislegend()
	fig
end

# ╔═╡ 0351835c-6de7-4230-ba53-c5c5921fbac9
p_total = 65.0

# ╔═╡ 40924d18-c8e8-4a6d-b58d-1797ce9b499a
data_mix = CSV.read("IRMOF-1_methane_ethane_mixture_isotherm_65bar_298K.csv", DataFrame)

# ╔═╡ ff136163-0e19-404f-9013-7f3f98f2d970
for gas in gases
	data_mix[:, "pred $gas [mmol/g]"] = [NaN for _ = 1:nrow(data_mix)]
end

# ╔═╡ 46024e58-2604-44d8-b413-afcc43d490fe
begin
	y = Dict("ethane" => collect(range(0.0, 0.12, length=25)[2:end]))
	y["methane"] = [1.0 - yᵢ for yᵢ in y["ethane"]]
	n = Dict(gas => zeros(24) for gas in gases)
end

# ╔═╡ e2a8e7f7-bcef-4dbf-ab85-60404c23c1ae
for i = 1:length(y["ethane"])
	partial_pressures = Dict(
		"ethane"  => p_total * y["ethane"][i], 
		"methane" => p_total * y["methane"][i]
	)
	a = iast(partial_pressures, models, p_warn=p_warn)
	for gas in gases
		n[gas][i] = a[gas]
	end
end

# ╔═╡ 9eba5791-71b2-46d7-880b-1badcf75ca1c
data_mix

# ╔═╡ 2abbcbe3-3dfb-4ee7-b8b0-864093088487
gas_to_color = Dict(zip(gases, ["red", "green"]))

# ╔═╡ 1fd0ed49-475a-4519-a15b-ca513bd4d7e1
function viz_mix(data_mix::DataFrame, y, n)
	fig = Figure()
	ax  = Axis(fig[1, 1], xlabel="y ethane", ylabel="mmol/g")
	scatter!(data_mix[:, "y_ethane"], data_mix[:, "EthaneLoading(mmol/g)"], 
		color=gas_to_color["ethane"])
	scatter!(data_mix[:, "y_ethane"], data_mix[:, "MethaneLoading(mmol/g)"], 
		color=gas_to_color["methane"])
	for gas in gases
		lines!(y["ethane"], n[gas], color=gas_to_color[gas])
	end
	return fig
end

# ╔═╡ e12b1eb0-9fc2-455d-836b-e7d14f7bc248
viz_mix(data_mix, y, n)

# ╔═╡ Cell order:
# ╠═8c4dd1f8-0fa8-11ed-2080-0fa9ccb116fa
# ╠═4e26b473-510c-4fbc-9bc6-283cf88e5f2a
# ╠═b1a58980-60da-4e75-b9b3-46ce0c31cab5
# ╠═ce6bcada-a06f-476b-8255-e651e60dfcec
# ╠═c43ae81a-8686-4237-b08b-06f3020cdd0b
# ╠═8b827180-af14-4590-8f23-13a219f84a72
# ╠═3b72d972-b4da-480f-beff-2dc1e670f04e
# ╠═ad050898-0d2f-4a6f-9c75-c4d093b16e04
# ╠═5b1cb128-4714-4024-bcc7-5fb4af9dff52
# ╠═0351835c-6de7-4230-ba53-c5c5921fbac9
# ╠═40924d18-c8e8-4a6d-b58d-1797ce9b499a
# ╠═ff136163-0e19-404f-9013-7f3f98f2d970
# ╠═46024e58-2604-44d8-b413-afcc43d490fe
# ╠═e2a8e7f7-bcef-4dbf-ab85-60404c23c1ae
# ╠═9eba5791-71b2-46d7-880b-1badcf75ca1c
# ╠═2abbcbe3-3dfb-4ee7-b8b0-864093088487
# ╠═1fd0ed49-475a-4519-a15b-ca513bd4d7e1
# ╠═e12b1eb0-9fc2-455d-836b-e7d14f7bc248
