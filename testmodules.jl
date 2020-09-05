include("md3simulation.jl")


using Main.MD3simulation.DVector
using Main.MD3simulation
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

w = World3D()
# populate the world with preset particles. Largest particle's radius is 2.0
test_populate_world!(w)
create_linkcell!(w)

@testset "world module" begin
    @testset "world population" begin
        @test w.Δ == 4.0
        @test w.nx == 5
        @test w.ny == 5
        @test w.nz == 4
        @test w.Lx == 20.0
        @test w.Ly == 20.0
        @test w.Lz == 16.0
        @test cellindex(w, w.p[1].r⃗⁰) == (3, 3, 2) # center
        @test cellindex(w, w.p[2].r⃗⁰) == (5, 3, 2) # radial top mid
        @test cellindex(w, w.p[3].r⃗⁰) == (1, 3, 2) # radial bot mid
        @test cellindex(w, w.p[4].r⃗⁰) == (3, 1, 2) # azimuthal left mid
        @test cellindex(w, w.p[5].r⃗⁰) == (3, 5, 2) # azimuthal right mid
        @test cellindex(w, w.p[6].r⃗⁰) == (5, 1, 2) # radial top / azimuthal left
        @test cellindex(w, w.p[7].r⃗⁰) == (5, 5, 2) # radial top / azimuthal right
        @test cellindex(w, w.p[8].r⃗⁰) == (1, 1, 2) # radial bot / azimuthal left
        @test cellindex(w, w.p[9].r⃗⁰) == (1, 5, 2) # radial bot / azimuthal right
        p = get_particle_image(w, w.p[1], 0, 0)
        @test (p.r⃗⁰.x == w.p[1].r⃗⁰.x) & (p.r⃗⁰.y == w.p[1].r⃗⁰.y) & (p.r⃗⁰.z == w.p[1].r⃗⁰.z)
        @test (p.r⃗¹.x == w.p[1].r⃗¹.x) & (p.r⃗¹.y == w.p[1].r⃗¹.y) & (p.r⃗¹.z == w.p[1].r⃗¹.z)
        p = get_particle_image(w, w.p[2], 1, 0)
        @test (p.r⃗⁰.x == w.p[2].r⃗⁰.x + w.Lx) & (p.r⃗⁰.y == w.p[2].r⃗⁰.y) & (p.r⃗⁰.z == w.p[2].r⃗⁰.z)
        @test (p.r⃗¹.x == w.p[2].r⃗¹.x) & (p.r⃗¹.y == w.p[2].r⃗¹.y-(3/2)*w.Ω₀*w.Lx) & (p.r⃗¹.z == w.p[2].r⃗¹.z)
        w.t = 2.0
        p = get_particle_image(w, w.p[3],-1, 0)
        @test (p.r⃗⁰.x == w.p[3].r⃗⁰.x - w.Lx) & (p.r⃗⁰.y == w.p[3].r⃗⁰.y + (3/2)*w.Lx*w.Ω₀*mod(w.t, w.Δt)) & (p.r⃗⁰.z == w.p[3].r⃗⁰.z)
        @test (p.r⃗¹.x == w.p[3].r⃗¹.x) & (p.r⃗¹.y == w.p[3].r⃗¹.y+(3/2)*w.Ω₀*w.Lx) & (p.r⃗¹.z == w.p[3].r⃗¹.z)
    end
end
