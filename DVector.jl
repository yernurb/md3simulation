module DVector

import Base: +, -, *, /, show
import LinearAlgebra: norm, ⋅, ×

abstract type Vec end

mutable struct Vec2<:Vec
    x::Float64
    y::Float64
end

mutable struct Vec3<:Vec
    x::Float64
    y::Float64
    z::Float64
end

vcopy(v::Vec2) = Vec2(v.x, v.y)
vcopy(v::Vec3) = Vec3(v.x, v.y, v.z)
norm(v::Vec2) = sqrt(v.x^2+v.y^2)
norm(v::Vec3) = sqrt(v.x^2+v.y^2+v.z^2)
norm2(v::Vec2) = v.x^2+v.y^2
norm2(v::Vec3) = v.x^2+v.y^2+v.z^2
+(v₁::Vec2, v₂::Vec2) = Vec2(v₁.x+v₂.x, v₁.y+v₂.y)
+(v₁::Vec3, v₂::Vec3) = Vec3(v₁.x+v₂.x, v₁.y+v₂.y, v₁.z+v₂.z)
-(v₁::Vec2, v₂::Vec2) = Vec2(v₁.x-v₂.x, v₁.y-v₂.y)
-(v₁::Vec3, v₂::Vec3) = Vec3(v₁.x-v₂.x, v₁.y-v₂.y, v₁.z-v₂.z)
*(a::Real, v::Vec2)  = Vec2(a*v.x, a*v.y)
*(a::Real, v::Vec3)  = Vec3(a*v.x, a*v.y, a*v.z)
*(v::Vec2, a::Real)  = Vec2(a*v.x, a*v.y)
*(v::Vec3, a::Real)  = Vec3(a*v.x, a*v.y, a*v.z)
/(v::Vec2, a::Real)  = Vec2(v.x/a, v.y/a)
/(v::Vec3, a::Real)  = Vec3(v.x/a, v.y/a, v.z/a)
⋅(v₁::Vec2, v₂::Vec2) = v₁.x*v₂.x + v₁.y*v₂.y
⋅(v₁::Vec3, v₂::Vec3) = v₁.x*v₂.x + v₁.y*v₂.y + v₁.z*v₂.z
×(v₁::Vec2, v₂::Vec2) = v₁.x*v₂.y-v₁.y*v₂.x
×(v₁::Vec3, v₂::Vec3) = Vec3(v₁.y*v₂.z-v₁.z*v₂.y, v₁.z*v₂.x-v₁.x*v₂.z, v₁.x*v₂.y-v₁.y*v₂.x)
setzero!(v::Vec2) = v.x=v.y=0
setzero!(v::Vec3) = v.x=v.y=v.z=0
unit(v::Vec2) = Vec2(v.x/norm(v), v.y/norm(v))
unit(v::Vec3) = Vec3(v.x/norm(v), v.y/norm(v), v.z/norm(v))

function add!(v::Vec2, u::Vec2)
    v.x += u.x
    v.y += u.y
end

function add!(v::Vec3, u::Vec3)
    v.x += u.x
    v.y += u.y
    v.z += u.z
end

function sub!(v::Vec2, u::Vec2)
    v.x -= u.x
    v.y -= u.y
end

function sub!(v::Vec3, u::Vec3)
    v.x -= u.x
    v.y -= u.y
    v.z -= u.z
end

function scale!(v::Vec2, a::Real)
    v.x *= a
    v.y *= a
end

function scale!(v::Vec3, a::Real)
    v.x *= a
    v.y *= a
    v.z *= a
end

function normalize!(v::Vec2)
    value = norm(v)
    v.x /= value
    v.y /= value
end

function normalize!(v::Vec3)
    value = norm(v)
    v.x /= value
    v.y /= value
    v.z /= value
end

function show(io::IO, v::Vec2)
    println(io, "($(v.x), $(v.y))")
end

function show(io::IO, v::Vec3)
    println(io, "($(v.x), $(v.y), $(v.z))")
end

export Vec2, Vec3, vcopy, norm, norm2, ⋅, ×, setzero!, unit, add!, sub!, scale!, normalize!

end # module
