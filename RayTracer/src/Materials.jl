module Materials

using Images
using FileIO
using LinearAlgebra

using ..GfxBase

export Material
export Lambertian, Specular, Dielectric, Glossy
export sample, pdf, brdfSample

mutable struct brdfSample
    omega::Vec3
    mult
end

abstract type Material end

mutable struct Lambertian <: Material
    albedo::RGB{Float32}
    specularColor::RGB{Float32} 
    specularExp
end

mutable struct Specular <: Material #Perfect mirror material
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

""" BRDF attenuation factor """
function eval(material::Lambertian, view::Vec3, omega::Vec3)
    return material.albedo * (1/pi)
end

function pdf(material::Lambertian, normal::Vec3, omega::Vec3)
    cosTheta = dot(normal, omega)
    #PDF for uniform hemisphere sampling is 1/2pi
    return cosTheta * (1/pi) #cosine weighted hemisphere pdf 
end

""" Sample the upper hemisphere around the surface normal. Assumes
    that view and normal are normalized """
function sample(material::Lambertian, view::Vec3, normal::Vec3)
    localOmega = cosineWeightedHemisphere()

    #Construct basis to bring omega into world space (the hemisphere aligned with the normal)
    w = normal

    u = normalize(cross(view,w))
    v = cross(w, u)

    omega = localOmega[1]*u + localOmega[2]*v + localOmega[3]*w
    cosTheta = dot(normal, omega)
    if (cosTheta < 0)
        omega = -omega
    end
    mult = (eval(material, view, omega) * cosTheta) * (1/pdf(material, normal, omega))
    return brdfSample(omega, mult)
end

""" Perfect specular material brdf """
function eval(Material::Specular, view::Vec3, omega::Vec3)
    return 0
end

function pdf(Material::Specular, normal::Vec3, omega::Vec3)
    return 0
end

""" Sample the upper hemisphere around the surface normal. Assumes
    that view and normal are normalized """
function sample(Material::Specular, view::Vec3, normal::Vec3)
    mirrorDirection = 2*dot(normal,view)*normal - view
    mult = RGB{Float32}(1.0, 1.0, 1.0)
    return brdfSample(mirrorDirection, mult)
end





end
