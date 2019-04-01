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
    @test quadgk(sin, 0, pi)[1]  ≈ 2
    @test abs(quadgk(x -> x^2, 0, 1)[2]) <= 1e-14

    ## SymPy
    @test limit(x -> sin(x)/x, 0)  ≈ 1
    @test integrate(sin,  0, pi)  ≈ 2



end

@testset "test functions" begin

    f(x) = sin(x)
    c = pi/4
    fn = tangent(f, c)
    @test fn(1)  ≈ f(c) + f'(c)*(1 - c)

    fn = secant(f, pi/6, pi/3)
    @test fn(pi/4) <= f(pi/4)

    out = lim(x -> sin(x)/x, 0)
    @test out[end, 2]  ≈ 1

    @test bisection(sin, 3, 4)  ≈ pi
    @test newton(sin, cos, 3.0)  ≈ pi
    @test newton(sin, 3.0)  ≈ pi

    @test riemann(sin, 0, pi, 10_000)  ≈ 2

end


@testset "2d" begin

    x = [[1,2,3], [4,5,6]]
    @test xs_ys(x)[1] == [1, 4]
    @test xs_ys(x)[2] == [2, 5]
    @test xs_ys(x)[3] == [3, 6]

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
        plot(xs_ys(r.(ts))...)
        arrow!(r(pi/4), r'(pi/4))

    end
end


end
