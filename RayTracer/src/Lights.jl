module Lights

using LinearAlgebra

export Light, DirectionalLight, AmbientLight, PointLight, AreaLight
export light_direction, shade_light
export sampleAreaLight, areaLightPDF, areaLightBRDF

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
    intensity::RGB{Float32}
    center::Vec3
    radius
end

""" Sample from the solid angle subtended by a spherical light source and return a point x' on the light's surface """
function sampleAreaLight(light::AreaLight, intersection::Vec3)
    c = light.center
    r = light.radius

    a = c - intersection #Direction from point towards light's center
    #Create an orthonormal basis where w is the vector from the point to the light's center, normalized
    w = normalize(a)
    #Create non-colinear vector t
    t = [w[1], w[2], w[3]]
    t[argmin(w)] = 1
    t = Vec3(t[1], t[2], t[3])
    
    u = normalize(cross(w,t))
    v = cross(w, u)

    #compute theta and phi values for sample in cone
    sinThetaMax2 = (r^2)/(dot(a, a))
    cosThetaMax = sqrt(max(0, 1 - sinThetaMax2))
    sp = Vec2(rand(), rand())
    cosTheta = (1 - sp[1]) + sp[1] * cosThetaMax
    sinTheta = sqrt(max(0, 1 - cosTheta^2))
    phi = 2*pi*sp[2]

    #Rotate point point psi from solid angle to cartesian coordinate system defined above
    psi = sinTheta*cos(phi)*u + sinTheta*sin(phi)*v + cosTheta*w
    #Find point x' on the light source
    lightRay = Ray(intersection, psi)
    xPrime = sphere_intersect(lightRay, light)
    return xPrime
end

""" Probability of having sampled x' """
function areaLightPDF(light::AreaLight, point::Vec3)
    sinThetaMax2 = (light.radius^2)/(dot(light.center - point, light.center - point))
    cosThetaMax = sqrt(max(0, 1 - sinThetaMax2))
    pdf = 1 / (2*pi*(1 - cosThetaMax))
    return pdf
end

""" Compute lighting via Blinn-Phong BRDF """
function areaLightBRDF(material, light::AreaLight, ray::Ray, normal::Vec3, intersection::Vec3, xp::Vec3)
    #Irradiance calculations
    light_vector = normalize(xp - intersection) #Treat sample point as if it were a singular point light
    distFromLight = norm(xp - intersection)
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
    
""" Calculate the direction of a given light source from the position point """
function light_direction(light::Light, point::Vec3) end

function light_direction(light::DirectionalLight, point::Vec3)
    return normalize(light.direction)
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

""" Ray-sphere intersection. """
function sphere_intersect(ray::Ray, light::AreaLight)
    squared_radius = light.radius * light.radius
    a = norm(ray.direction)*norm(ray.direction)
    half_b = dot(ray.direction, (ray.origin - light.center))
    c = -squared_radius + norm(ray.origin - light.center)*norm(ray.origin - light.center)
    discriminant = (half_b*half_b) - c*a
    if discriminant < 0
        return nothing
    else
        t_0 = (-half_b - sqrt(discriminant))/a
        if t_0 > 0
            intersection_pt = ray.origin + t_0 * ray.direction
        else
            t_1 = (-half_b + sqrt(discriminant))/a
            intersection_pt = ray.origin + t_1 * ray.direction
        end
        return intersection_pt
    end
end

end # module Lights
