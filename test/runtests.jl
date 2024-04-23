using MTH229


using Test


## test package
@testset "test packages" begin
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

    f₁ = x ->  1 + 100x^2 - x^3
    @test isempty(find_zeros(f₁, -100, 100) )
    f₁s = find_zeros(fisheye(f₁), -pi/2, pi/2)
    @test !isempty(f₁s)
    @test abs(f₁(tan(first(f₁s)))) <= 1e-4


    @test bisection(sin, 3, 4)  ≈ pi
    @test newton(sin, cos, 3.0)  ≈ pi
    @test newton(sin, 3.0)  ≈ pi

    @test fubini((x,y) -> 1, (x->-sqrt(1-x^2), x->sqrt(1-x^2)), (-1,1)) ≈ pi

end
