push!(LOAD_PATH, joinpath(pwd(), "../src"))

using IAST, Test, LinearAlgebra

@testset "binary Langmuir tests" begin
    a_true(p, i, K, M) = M * K[i] * p[i] / (1 + dot(K, p))

    M = 3.0
    K = [2.2, 5.6]
    aims = [LangmuirModel(M, K[i]) for i = 1:2]
    
    for i = 1:25
        p = rand(2)

        a = iast(p, aims)
        @test all([a_true(p, i, K, M) for i = 1:2] .≈ a)
    end
end
