push!(LOAD_PATH, joinpath(pwd(), "../src"))

using IAST, Test, LinearAlgebra, DataFrames

@testset "binary Langmuir tests" begin
    a_true(p, i, K, M) = M * K[i] * p[i] / (1 + dot(K, p))

    M = 3.0
    K = [2.2, 5.6]
    aims = [LangmuirModel(M, K[i]) for i = 1:2]
    
    for i = 1:100
        p = rand(2)

        a = iast(p, aims)
        @test all([a_true(p, i, K, M) for i = 1:2] .≈ a)
    end
end

@testset "isotherm fit tests" begin
    LangmuirModel()
	my_data = DataFrame(p=range(0.0, 100.0, length=10000))
	my_data[:, "n"] = 3 * my_data[:, "p"] ./ (1 .+ my_data[:, "p"])

	θ₀ = IAST._default_θ_guess(LangmuirModel(), my_data, "p", "n")

	res = fit(my_data, "p", "n", LangmuirModel(M=0.1, K=2.0))

	@test isapprox(θ₀.K, 1.0, atol=0.01)
	@test isapprox(θ₀.M, 3.0, atol=0.05)
	@test isapprox(res.K, 1.0, atol=0.01)
	@test isapprox(res.M, 3.0, atol=0.01)
end
