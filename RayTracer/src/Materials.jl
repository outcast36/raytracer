module Materials

using Images
using FileIO

using ..GfxBase
using ..Lights

export Material
export Lambertian, Metallic, Dielectric, Glossy

abstract type Material end

mutable struct Lambertian <: Material
    albedo::RGB{Float32}
    specularColor::RGB{Float32} 
    specularExp
end

mutable struct Metallic <: Material
    albedo::RGB{Float32}
    specularColor::RGB{Float32}
    specularExp
end

mutable struct Dielectric <: Material
    index::Float64
end

mutable struct Glossy <: Material
    albedo::RGB{Float32}
    specularColor::RGB{Float32} 
    reflectance::Float64 
    specularExp
end

#TODO: Make translucent material to keep ideal dielectric sepearate from glossy transparency?

end
