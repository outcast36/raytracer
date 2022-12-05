module Materials

using Images
using FileIO

using ..GfxBase
using ..Lights

export Material
export Lambertian, Metallic, Dielectric

abstract type Material end

mutable struct Lambertian <: Material
    albedo::RGB{Float32}
    specularColor::RGB{Float32}
    specularExp
end

mutable struct Metallic <: Material
    albedo::RGB{Float32}
    specularColor::RGB{Float32}
    reflectance::Float64
    specularExp
end

mutable struct Dielectric <: Material
    index::Float64
end

end