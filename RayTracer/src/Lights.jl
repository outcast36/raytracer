module Lights

using LinearAlgebra

export Light, DirectionalLight, AmbientLight, PointLight
export light_direction, shade_light

using ..GfxBase

# Light types:
#  - directional is a distant source whose direction is always the same
#  - point light emits light from a single position in 3D
abstract type Light end
struct DirectionalLight <: Light
    normalIrradiance::RGB{Float32}
    direction::Vec3
end

struct PointLight <: Light
    intensity::RGB{Float32}
    position::Vec3
end

""" Model area light source as a sphere """
struct AreaLight <: Light
    center::Vec3
    intensity::RGB{Float32}
    radius
end

""" Sample a spherical light-source on the hemisphere visible from the point """
function sampleAreaLight(light::AreaLight, point::Vec3)
    radius = light.radius
    sp = Vec2(rand(), rand()) #Two canonical random numbers
    if (norm(light.center - point)^2 < radius^2)
        return uniformSampleSphere(sp) #point inside light source
    else
        #Create a basis where w is the vector from the light's center to the point
        w = normalize(light.center - point)
        #Create non-colinear vector t
        t = [w[1], w[2], w[3]]
        t[argmin(w)] = 1
        t = Vec3(t[1], t[2], t[3])
        
        u = normalize(cross(t,w))
        v = cross(w, u)

        #compute theta and phi values for sample in cone
        sinThetaMax2 = (radius^2)/(dot(light.center - point, light.center - point))
        cosThetaMax = sqrt(max(0, 1 - sinThetaMax2))
        cosTheta = (1 - sample[1]) + sp[1] * cosThetaMax
        sinTheta = sqrt(max(0, 1 - cosTheta^2))
        phi = 2*pi*sp[2]

        #distTocenter = norm(light.center - point)
        #distToSurface = distTocenter * cosTheta - sqrt(radius^2 - (distTocenter^2 * sinTheta^2))
        #cosAlpha = (distTocenter^2+radius^2-distToSurface^2)/(2*radius*distTocenter)
        #sinAlpha = sqrt(max(0, 1 - cosAlpha^2))
        #compute angle alpha from center of sphere to sampled point on surface
        normal = -sinTheta*cos(phi)*u - sinTheta*sin(phi)*v - cosTheta*w
        sampledLight = radius * normal
        #compute surface normal and sampled point on sphere
        return sampledLight
    end
end

""" Calculate the direction of a given light source from the position point """
function light_direction(light::Light, point::Vec3) end

function light_direction(light::DirectionalLight, point::Vec3)
    return normalize(light.direction)
end

function light_direction(light::AreaLight, point::Vec3)
    samplePoint = sampleAreaLight(light, point)#randomly sample point on light geometry
    return normalize(samplePoint - point) #direction from point to light source; normalized
end

# point source: direction from the point to the light position, normalized
function light_direction(light::PointLight, point::Vec3)
    return normalize(light.position - point)
end

""" Determine point-light contribution for shading a diffuse surface """
function shade_light(material, ray::Ray, normal::Vec3, intersection::Vec3, light::PointLight)
    #Irradiance calculations
    light_vector = light_direction(light, intersection)
    distFromLight = norm(intersection - light.position)
    n_dot_l = dot(normal, light_vector)
    irradiance = light.intensity * (1/(distFromLight^2)) * max(0, n_dot_l) #RGB color

    #!! -- RGB color * RGB color -- !!#
    diffuse_light = colorMultiply(material.albedo, irradiance)

    #Specular light calculations
    view_direction = -ray.direction
    half_vector = normalize(light_vector + view_direction)
    n_dot_h = dot(normal, half_vector)

    #!! -- RGB color * RGB color -- !!#
    specular_light = (max(0, n_dot_h)^material.specularExp) * colorMultiply(material.specularColor, irradiance)
    
    return diffuse_light + specular_light #RGB color + RGB color
end

""" Determine direction-light contribution for shading a diffuse surface """
function shade_light(material, ray::Ray, normal::Vec3, intersection::Vec3, light::DirectionalLight)
    #Irradiance calculations
    light_vector = light_direction(light, intersection)
    n_dot_l = dot(normal, light_vector)
    irradiance = light.normalIrradiance * max(0, n_dot_l)  #RGB color
    
    #!! -- RGB color * RGB color -- !!#
    diffuse_light = colorMultiply(material.albedo, irradiance)
    
    #Specular light calculations
    view_direction = -ray.direction
    half_vector = normalize(light_vector + view_direction)
    n_dot_h = dot(normal, half_vector)

    #!! -- RGB color * RGB color -- !!#
    specular_light = (max(0, n_dot_h)^material.specularExp) * colorMultiply(material.specularColor, irradiance)
    
    return diffuse_light + specular_light #RGB color + RGB color
end

""" Determine area light contribution for shading a diffuse surface by sampling random points on the spherical light source """
function shade_light(material, ray::Ray, normal::Vec3, intersection::Vec3, light::AreaLight)
    #Irradiance calculations
    light_vector = light_direction(light, intersection)
    distFromLight = norm(intersection - light.position)
    n_dot_l = dot(normal, light_vector)
    irradiance = light.intensity * (1/(distFromLight^2)) * max(0, n_dot_l) #RGB color

    #!! -- RGB color * RGB color -- !!#
    diffuse_light = colorMultiply(material.albedo, irradiance)

    #Specular light calculations
    view_direction = -ray.direction
    half_vector = normalize(light_vector + view_direction)
    n_dot_h = dot(normal, half_vector)

    #!! -- RGB color * RGB color -- !!#
    specular_light = (max(0, n_dot_h)^material.specularExp) * colorMultiply(material.specularColor, irradiance)
    
    return diffuse_light + specular_light #RGB color + RGB color
end

end # module Lights
