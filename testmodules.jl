include("DVector.jl")

using Main.DVector
using Test

vec2_1 = Vec2( 2.0, 3.0)
vec2_2 = Vec2( 1.0, 2.0)
vec2_3 = Vec2(-1.0, 4.5)

vec3_1 = Vec3( 2.0, 3.0, 0.0)
vec3_2 = Vec3( 1.0, 2.0, 1.5)
vec3_3 = Vec3(-1.0, 4.5,-2.5)

@testset "custom vector module" begin
    @testset "2D vectors" begin
        v_test = vec2_1 + vec2_2
        @test (v_test.x == 3.0) & (v_test.y == 5.0)
        v_test = vec2_1 + vec2_2 + vec2_3
        @test (v_test.x == 2.0) & (v_test.y == 9.5)
        v_test = vec2_1 - vec2_2
        @test (v_test.x == 1.0) & (v_test.y == 1.0)
        v_test = vec2_1 - vec2_2 - vec2_3
        @test (v_test.x == 2.0) & (v_test.y == -3.5)
        @test norm(vec2_1) ≈ √13
        v_test = 3*vec2_3
        @test (v_test.x == -3.0) & (v_test.y == 13.5)
        v_test = vec2_3*3
        @test (v_test.x == -3.0) & (v_test.y == 13.5)
        v_test = vec2_3 / 2
        @test (v_test.x == -0.5) & (v_test.y == 2.25)
        @test vec2_1⋅vec2_2 == 8.0
        @test vec2_1⋅vec2_3 == 11.5
        @test vec2_1×vec2_2 == 1.0
        @test vec2_1×vec2_3 == 12.0
        v_test = vcopy(vec2_2)
        setzero!(v_test)
        @test (v_test.x == 0) & (v_test.y == 0)
        v_test = unit(vec2_2)
        @test (v_test.x ≈ 1.0/√5) & (v_test.y ≈ 2.0/√5)
        v_test = vcopy(vec2_1)
        add!(v_test, vec2_2)
        @test (v_test.x == 3.0) & (v_test.y == 5.0)
        v_test = vcopy(vec2_1)
        sub!(v_test, vec2_2)
        @test (v_test.x == 1.0) & (v_test.y == 1.0)
        v_test = vcopy(vec2_3)
        scale!(v_test, 3)
        @test (v_test.x == -3.0) & (v_test.y == 13.5)
        v_test = vcopy(vec2_2)
        normalize!(v_test)
        @test (v_test.x ≈ 1.0/√5) & (v_test.y ≈ 2.0/√5)
    end
    @testset "3D vectors" begin
        v_test = vec3_1 + vec3_2
        @test (v_test.x == 3.0) & (v_test.y == 5.0) & (v_test.z == 1.5)
        v_test = vec3_1 + vec3_2 + vec3_3
        @test (v_test.x == 2.0) & (v_test.y == 9.5) & (v_test.z == -1.0)
        v_test = vec3_1 - vec3_2
        @test (v_test.x == 1.0) & (v_test.y == 1.0) & (v_test.z == -1.5)
        v_test = vec3_1 - vec3_2 - vec3_3
        @test (v_test.x == 2.0) & (v_test.y == -3.5) & (v_test.z == 1.0)
        @test norm(vec3_3) ≈ √27.5
        v_test = 3*vec3_3
        @test (v_test.x == -3.0) & (v_test.y == 13.5) & (v_test.z == -7.5)
        v_test = vec3_3*3
        @test (v_test.x == -3.0) & (v_test.y == 13.5) & (v_test.z == -7.5)
        v_test = vec3_3/2
        @test (v_test.x == -0.5) & (v_test.y == 2.25) & (v_test.z == -1.25)
        @test vec3_1⋅vec3_2 == 8.0
        @test vec3_2⋅vec3_3 == 4.25
        v_test = vec3_1×vec3_2
        @test (v_test.x == 4.5) & (v_test.y == -3.0) & (v_test.z == 1.0)
        v_test = vec3_2×vec3_3
        @test (v_test.x == -11.75) & (v_test.y == 1.0) & (v_test.z == 6.5)
        v_test = vcopy(vec3_2)
        setzero!(v_test)
        @test (v_test.x == 0) & (v_test.y == 0) & (v_test.z == 0)
        v_test = unit(vec3_2)
        @test (v_test.x == 1.0/√7.25) & (v_test.y == 2.0/√7.25) & (v_test.z == 1.5/√7.25)
        v_test = vcopy(vec3_1)
        add!(v_test, vec3_2)
        @test (v_test.x == 3.0) & (v_test.y == 5.0) & (v_test.z == 1.5)
        v_test = vcopy(vec3_1)
        sub!(v_test, vec3_2)
        @test (v_test.x == 1.0) & (v_test.y == 1.0) & (v_test.z == -1.5)
        v_test = vcopy(vec3_3)
        scale!(v_test, 3)
        @test (v_test.x == -3.0) & (v_test.y == 13.5) & (v_test.z == -7.5)
        v_test = vcopy(vec3_2)
        normalize!(v_test)
        @test (v_test.x == 1.0/√7.25) & (v_test.y == 2.0/√7.25) & (v_test.z == 1.5/√7.25)
    end
end