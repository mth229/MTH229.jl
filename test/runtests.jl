using MTH229


using Test


## test package
@testset "test packages" begin

    ## Roots
    @test fzero(sin, 3, 4)  ≈ pi
    @test fzero(sin, 3.0)  ≈ pi

    ## ForwardDiff
    f(x) = sin(x)
    @test f'(2)  ≈ cos(2)
    @test f''(2)  ≈ -sin(2)

    ## QuadGK
    @test first(quadgk(sin, 0, pi))  ≈ 2
    @test abs(quadgk(x -> x^2, 0, 1)[2]) <= 1e-14

    ## SymPy
    @syms x
    @test limit(sin(x)/x, x=>0) == 1
    @test integrate(sin(x),  (x, 0, pi)) == 2

    ## SimpleExpressions
    @symbolic x
    @test solve(sin(x) ~ 0, (3,4)) ≈ pi
    @test solve(sin(x) ~ 0, 3) ≈ pi
    @test solve(sin(x) ~ 0, 3..4) ≈ [pi]


end

@testset "test functions" begin

    f(x) = sin(x)
    c = pi/4
    fn = tangent(f, c)
    @test fn(1)  ≈ f(c) + f'(c)*(1 - c)

    fn = secant(f, pi/6, pi/3)
    @test fn(pi/4) <= f(pi/4)

    out = lim(x -> sin(x)/x, 0)
    @test_broken out[end, 2]  ≈ 1 # an iterator now

    @test bisection(sin, 3, 4)  ≈ pi
    @test newton(sin, cos, 3.0)  ≈ pi
    @test newton(sin, 3.0)  ≈ pi

    f₁ = x ->  1 + 100x^2 - x^3
    @test isempty(find_zeros(f₁, -100, 100) )
    #f₁s = find_zeros(fisheye(f₁), -pi/2, pi/2) # need to export
    #@test !isempty(f₁s)
    #@test abs(f₁(tan(first(f₁s)))) <= 1e-4

    out = sign_chart(x -> (x-1)*(x-2)/(x-3), 0, 4)
    @test all([o[1] for o ∈ out] .≈[1,2,3])

    @test riemann(sin, 0, pi, 10_000)  ≈ 2
    @test_broken fubini((x,y) -> 1, (x->-sqrt(1-x^2), x->sqrt(1-x^2)), (-1,1)) ≈ pi # need to export
end


@testset "2d" begin

    x = [[1,2,3], [4,5,6]]
    @test unzip(x)[1] == [1, 4]
    @test unzip(x)[2] == [2, 5]
    @test unzip(x)[3] == [3, 6]

    @test length(unzip(x -> x, 0, 1)[1])  <= 50 # 21
    @test length(unzip(x-> sin(10pi*x), 0, 1)[1]) >= 50 # 233

    @test uvec([2,2]) == 1/sqrt(2) * [1,1]

end

@testset "plots" begin

    if isinteractive()

        f(x) = 1/x

        trimplot(f, -1, 1)

        plotif(f, f, -1, 1)
        plotif(f, f', -1, 1)

        signchart(sin, 0, 4pi)

        newton_vis(log, 1/2)

        r(t) = [sin(t), cos(t)]
        ts = range(0, stop=pi/2, length=100)
        plot(unzip(r.(ts))...)
        arrow!(r(pi/4), r'(pi/4))

        V(x,y) = [x, x-y]
        plot()
        vectorfieldplot!(V)

    end

end
