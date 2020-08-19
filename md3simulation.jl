module MD3simulation

include("DVector.jl")
using .DVector

struct Particle
    r⃗⁰::Vec3   # position vector (zero's order time derivative)
    r⃗¹::Vec3   # velocity vector (first order time derivative)
    r⃗²::Vec3   # acceleartion vector (second order time derivative)
    r⃗³::Vec3   # third order time derivative of position vector
    r⃗⁴::Vec3   # fourth order time derivative of position vector
    F⃗::Vec3    # actual resulting force acting upon the particle
    M::Float64 # mass
    R::Float64 # radius
    κ::Float64 # elastic coefficient
    γ::Float64 # damping coefficient
    cc::Set{Int} # current contacts
    function Particle(mass = 0.0, radius = 0.0, elastic = 0.0, viscous = 0.0)
        @assert(0 ≤ mass, "Particle mass should be non-negative")
        @assert(0 ≤ radius, "Particle radius should be non-negative")
        @assert(0 ≤ elastic, "Elasticity coeffitient should be non-negative")
        @assert(0 ≤ viscous, "Viscous damping coefficient should be non-negative")
        r⃗⁰ = Vec3(0, 0, 0)
        r⃗¹ = Vec3(0, 0, 0)
        r⃗² = Vec3(0, 0, 0)
        r⃗³ = Vec3(0, 0, 0)
        r⃗⁴ = Vec3(0, 0, 0)
        F⃗ = Vec3(0, 0, 0)
        cc = Set{Int}()
        new(r⃗⁰, r⃗¹, r⃗², r⃗³, r⃗⁴, F⃗, mass, radius, elastic, viscous, cc)
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
    nb::Vector{Vector{Int}} # collection of neighbouring indices
    sd::Vector{Any} # system data of various types, to test the system
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
        new(Lx,Ly,Lz,Ω₀,t,δt,Δt,p,Δ,nx,ny,nz,lc,nb,sd)
    end
end

# TODO standard population method

# populate the world with preset particles at preset position for further testing
function test_populate_world!(w::World3D)
    w.Lx = 20.0
    w.Ly = 18.0
    w.Lz = 15.0
    w.Ω₀ = 2.0
    w.Δt = (2/3w.Ω₀)*(w.Ly/w.Lx)
    # particle 1 at the center
    p = Particle(1.0, 1.0, 1.0, 1.0)
    p.r⃗⁰.x = w.Lx/2; p.r⃗⁰.y = w.Ly/2; p.r⃗⁰.z = w.Lz/2
    p.r⃗¹.x = 1.0; p.r⃗¹.y =-1.0; p.r⃗¹.z = 0.0
    push!(w.p, p)
    # particle 2 at the radial top edge
    p = Particle(1.0, 2.0, 1.0, 1.0)
    p.r⃗⁰.x = w.Lx-2.0; p.r⃗⁰.y = w.Ly/2; p.r⃗⁰.z = w.Lz/2
    p.r⃗¹.x =-1.0; p.r⃗¹.y = 2.1; p.r⃗¹.z = 0.2
    push!(w.p, p)
    # particle 3 at the radial bottom edge
    p = Particle(1.0, 1.5, 1.0, 1.0)
    p.r⃗⁰.x = 1.5; p.r⃗⁰.y = w.Ly/2; p.r⃗⁰.z = w.Lz/2
    p.r⃗¹.x = 1.0; p.r⃗¹.y = 0.1; p.r⃗¹.z =-0.2
    push!(w.p, p)
    # particle 4 at the azimuthal left edge
    p = Particle(1.0, 1.0, 1.0, 1.0)
    p.r⃗⁰.x = w.Lx/2; p.r⃗⁰.y = 1.0; p.r⃗⁰.z = w.Lz/2
    p.r⃗¹.x =-1.0; p.r⃗¹.y = 0.5; p.r⃗¹.z = 0.2
    push!(w.p, p)
    # particle 5 at the azimuthal right edge
    p = Particle(1.0, 2.0, 1.0, 1.0)
    p.r⃗⁰.x = w.Lx/2; p.r⃗⁰.y = w.Ly-2.0; p.r⃗⁰.z = w.Lz/2
    p.r⃗¹.x = 1.5; p.r⃗¹.y =-1.2; p.r⃗¹.z = 2.2
    push!(w.p, p)
    # particle 6 at the radial top and azimuthal left corner
    p = Particle(1.0, 1.0, 1.0, 1.0)
    p.r⃗⁰.x = w.Lx-1.0; p.r⃗⁰.y = 1.0; p.r⃗⁰.z = w.Lz/2
    p.r⃗¹.x =-0.8; p.r⃗¹.y = 1.1; p.r⃗¹.z = 0.2
    push!(w.p, p)
    # particle 7 at the radial top and azimuthal right corner
    p = Particle(1.0, 1.4, 1.0, 1.0)
    p.r⃗⁰.x = w.Lx-1.4; p.r⃗⁰.y = w.Ly-1.0; p.r⃗⁰.z = w.Lz/2
    p.r⃗¹.x = 1.8; p.r⃗¹.y = 3.1; p.r⃗¹.z =-0.2
    push!(w.p, p)
    # particle 8 at the radial bottom and azimuthal left corner
    p = Particle(1.0, 0.2, 1.0, 1.0)
    p.r⃗⁰.x = 0.2; p.r⃗⁰.y = 0.2; p.r⃗⁰.z = w.Lz/2
    p.r⃗¹.x =-2.1; p.r⃗¹.y =-1.1; p.r⃗¹.z =-0.2
    push!(w.p, p)
    # particle 9 at the radial bottom and azimuthal right corner
    p = Particle(1.0, 1.2, 1.0, 1.0)
    p.r⃗⁰.x = 1.2; p.r⃗⁰.y = w.Ly-1.2; p.r⃗⁰.z = w.Lz/2
    p.r⃗¹.x = 2.0; p.r⃗¹.y = 0.1; p.r⃗¹.z = 0.2
    push!(w.p, p)
end


"""
    create_linkcell!(w::World3D)

Creates the linkcell structure based on the particles' sizes in the world.
Sets the size of the cell equal to the diameter of the largest particle in the world.
"""
function create_linkcell!(w::World3D)
    Δ = 2maximum(collect(p.R for p in w.p)) # cell size is equal to the diameter of the largest particle
    w.Δ = Δ
    w.nx = ceil(w.Lx/Δ) # estimate number of cells in x direction
    w.ny = ceil(w.Ly/Δ) # estimate number of cells in y direction
    w.nz = ceil(w.Lz/Δ) # estimate number of cells in z direction
    w.Lx = w.nx*Δ      # adjust the sizes of the world to better tiling
    w.Ly = w.ny*Δ
    w.Lz = w.nz*Δ
    w.Δt = (2/3w.Ω₀)*(w.Ly/w.Lx)

    w.lc = Array{Set{Int64}}(undef, w.nx, w.ny, w.nz) # create linkcell structure as a 3D array of sets
    for i in CartesianIndices(w.lc)          # this structure will hold the
        w.lc[i] = Set{Int64}()               # particle id's in corresponding cells
    end

    for i in eachindex(w.p)
        cx, cy, cz = cellindex(w, w.p[i].r⃗⁰)
        push!(w.lc[cx, cy, cz], i)
    end
    nothing
end

# get cell index of a given vector
function cellindex(w::World3D, r⃗::Vec3)
    return (Int(w.nx*r⃗.x÷w.Lx)+1, Int(w.ny*r⃗.y÷w.Ly)+1, Int(w.nz*r⃗.z÷w.Lz)+1)
end

# TODO find images of a particle in XY plane
function get_particle_image(w::World3D, p::Particle, n::Int, m::Int)
    τ = mod(w.t, w.Δt)
    p_image = Particle()
end

# TODO find neighbouring cells, taking the sheared borders into account
# TODO check distance between two particles and their images
# TODO calculate forces between particles
# TODO predictor
# TODO corrector
# TODO make step and update world


export Particle, World3D, test_populate_world!, create_linkcell!, cellindex

end  # module MD3simulation
