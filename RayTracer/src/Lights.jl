module Lights

using LinearAlgebra

export Light, SphericalLight
export sampleLight, lightPDF

using ..GfxBase

abstract type Light end

""" Model area light source as a sphere """
struct SphericalLight <: Light
    geometry
end

""" Model area light with a single triangle """
struct TriangleLight <: Light
    geometry
end

#disk light

#cylinder/surface of rotation light

""" Sample a point on a spherical light in the solid angle subtended by the light at point x """
function sampleLight(light::SphericalLight, intersection::Vec3)
    #PHYSICALLY CORRECT DIRECT LIGHTING FOR DISTRIBUTION RAY TRACING, Page 4 -- Changyaw Wang, 1992
    c = light.geometry.center
    r = light.geometry.radius

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

""" Sample a point uniformly over a given triangle's area """
function sampleLight(light::TriangleLight, x::Vec3)
    return 0
end

""" Probability of having sampled x' with solid angle sampling """
function lightPDF(light::SphericalLight, point::Vec3, xp::Vec3)
    c = light.geometry.center
    r = light.geometry.radius

    lightNormal = (1/r) * (xp - c) #Direction from light center to sampled point x'
    omega = normalize(xp-point) #Direction from intersection point to the point x' on the luminaire
    cosThetaPrime = max(0, dot(lightNormal, -omega))       
    distFactor = 1/(norm(xp - point)*norm(xp-point))
    
    sinThetaMax2 = (r*r)/(norm(c - point)*norm(c - point))
    cosThetaMax = sqrt(max(0, 1 - sinThetaMax2))
    pdf = (1 / (2*pi*(1 - cosThetaMax))) * cosThetaPrime*distFactor
    return pdf
end

""" Ray-sphere intersection. """
function sphere_intersect(ray::Ray, light::SphericalLight)
    lightCenter = light.geometry.center
    squared_radius = light.geometry.radius * light.geometry.radius
    a = norm(ray.direction)*norm(ray.direction)
    half_b = dot(ray.direction, (ray.origin - lightCenter))
    c = -squared_radius + norm(ray.origin - lightCenter)*norm(ray.origin - lightCenter)
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