push!(LOAD_PATH, joinpath(pwd(), "../src"))

using IAST, Test, LinearAlgebra, DataFrames, CSV

@testset "compare to mixed-gas sim tests" begin
    data = Dict(gas => CSV.read("IRMOF-1_$(gas)_isotherm_298K.csv", DataFrame) 
                    for gas in ["methane", "ethane"])
    data_mix = CSV.read("IRMOF-1_methane_ethane_mixture_isotherm_65bar_298K.csv", DataFrame)
end

@testset "binary Langmuir tests" begin
    a_true(p, i, K, M) = M * K[i] * p[i] / (1 + dot(K, p))

    M = 3.0
    K = [2.2, 5.6]
    aims = [LangmuirModel(M, K[i]) for i = 1:2]
    
    for i = 1:100
        p = 0.0001 .+ rand(2)

        a = iast(p, aims)
        @test all([a_true(p, i, K, M) for i = 1:2] .≈ a)
    end
end

@testset "isotherm tests" begin
    K = 1 + rand()
    M = 1 + rand()

    lm = LangmuirModel(K=K, M=M)
    tm = TemkinApproxModel(K=K, M=M, θ=0.0)
    qm = QuadraticModel(K=K, M=M/2, ϕ=1.0)
    lm2 = LangmuirModel(K=K + rand(), M=M + rand())
    dslm = DualSiteLangmuirModel(M₁=lm.M, K₁=lm.K, M₂=lm2.M, K₂=lm2.K)

    p = 1 + rand()

    @test isapprox(loading(p, lm), loading(p, tm))
    @test isapprox(loading(p, lm) + loading(p, lm2), loading(p, dslm))
    @test isapprox(loading(p, lm), loading(p, qm))
    @test isapprox(grand_pot(p, lm), grand_pot(p, tm))
    @test isapprox(grand_pot(p, lm), grand_pot(p, qm))
    @test isapprox(grand_pot(p, lm) + grand_pot(p, lm2), grand_pot(p, dslm))
end

@testset "isotherm fit tests" begin
    LangmuirModel()
	my_data = DataFrame(p=range(0.0, 100.0, length=10000))
	my_data[:, "n"] = 3 * my_data[:, "p"] ./ (1 .+ my_data[:, "p"])
    ads_data = AdsorptionIsothermData(my_data, "p", "n")

    θ₀ = IAST._default_θ_guess(ads_data, LangmuirModel())

	res = identify_params(ads_data, LangmuirModel(M=0.1, K=2.0))

	@test isapprox(θ₀.K, 1.0, atol=0.01)
	@test isapprox(θ₀.M, 3.0, atol=0.05)
	@test isapprox(res.K, 1.0, atol=0.01)
	@test isapprox(res.M, 3.0, atol=0.01)
end

