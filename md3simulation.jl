module MD3simulation

include("DVector.jl")

using Main.DVector

struct Particle
    r⃗⁰::Vec3   # position vector
    r⃗¹::Vec3   # velocity vector
    r⃗²::Vec3   # acceleartion vector
    r⃗³::Vec3   # third time derivative of position vector
    r⃗⁴::Vec3   # fourth time derivative of position vector
    F⃗::Vec3    # actual force acting upon the particle
    M::Float64 # mass
    R::Float64 # radius
    Y::Float64 # elastic coefficient
    γ::Float64 # damping coefficient
    cc::Set{Int} # current contacts
    function Particle(mass = 0.0, radius = 0.0, young = 0.0, damping = 0.0)
        @assert(0 ≤ mass, "Particle mass should be non-negative")
        @assert(0 ≤ radius, "Particle radius should be non-negative")
        @assert(0 ≤ young, "Young modulus should be non-negative")
        @assert(0 ≤ damping, "Damping coefficient should be non-negative")
        r⃗⁰ = Vec3(0, 0, 0)
        r⃗¹ = Vec3(0, 0, 0)
        r⃗² = Vec3(0, 0, 0)
        r⃗³ = Vec3(0, 0, 0)
        r⃗⁴ = Vec3(0, 0, 0)
        F⃗ = Vec3(0, 0, 0)
        cc = Set{Int}()
        new(r⃗⁰, r⃗¹, r⃗², r⃗³, r⃗⁴, F⃗, mass, radius, young, damping, cc)
    end
end


mutable struct World3D
    Lx::Float64 # radial size of the simulation area
    Ly::Float64 # azimuthal size of the simulation area
    Lz::Float64 # vertical size of the simulation area
    Ω₀::Float64 # orbital speed of the Hill's box
    t::Float64  # global simulation time
    δt::Float64 # time step of integration scheme
    Δt::Float64 # period of shearing motion
    p::Vector{Particle} # list of all particles
    Δ::Float64 # cell size of a cubic cell
    nx::Int64 # number of cells in x direction
    ny::Int64 # number of cells in y direction
    nz::Int64 # number of cells in z direction
    lc::Array{Set{Int}, 3} # linkcell structure
    sd::Vector{Any} # system data of various types, to test the system
    nb::Vector{Vector{Int}} # collection of neighbouring indices
    function World3D(rad_dim=1000, azi_dim=1000, vert_dim=1000, omega=1e-4)
        Lx = rad_dim
        Ly = azi_dim
        Lz = vert_dim
        Ω₀ = omega
        t = 0.0
        δt = 0.001
        p = Vector{Particle}()
        Δ = 1.0
        Δt = (2/3Ω₀)*(Ly/Lx)
        nx = Int(Lx÷Δ)
        ny = Int(Ly÷Δ)
        nz = Int(Lz÷Δ)
        Lx = nx*Δ
        Ly = ny*Δ
        Lz = nz*Δ
        lc = Array{Set{Int}, 3}(undef, nx, ny, nz)
        sd = Vector{Any}()
        nb = [[ 0, 0, 0],[-1, 1, 0],[ 0, 1, 0],[ 1, 1, 0],[ 1, 0, 0],
              [ 0, 0,-1],[-1, 1,-1],[ 0, 1,-1],[ 1, 1,-1],[ 1, 0,-1],
              [ 0, 0, 1],[-1, 1, 1],[ 0, 1, 1],[ 1, 1, 1],[ 1, 0, 1]]
        new(Lx,Ly,Lz,Ω₀,t,δt,Δt,p,Δ,nx,ny,nz,lc,sd,nb)
    end
end

# TODO standard population method
# TODO create linkcell structure for given population
# TODO get particle cell index
# TODO find images of a particle in XY plane
# TODO find neighbouring cells, taking the sheared borders into account
# TODO check distance between two particles and their images
# TODO calculate forces between particles
# TODO predictor
# TODO corrector
# TODO make step and update world

end  # module MD3simulation
