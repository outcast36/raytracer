""" 
Wyatt Ayers
November 23, 2022
Main module for raytracer
"""

module RayTracer

export main

using FileIO
using Images
using StaticArrays
using LinearAlgebra

push!(LOAD_PATH, pwd())
include("GfxBase.jl")
include("Lights.jl")
include("Materials.jl")
include("MeshGen.jl")
include("Scenes.jl")
include("Cameras.jl")
include("TestScenes.jl")

using .GfxBase
using .Lights
using .Materials

import .Scenes
import .Scenes.Scene
import .Scenes.HitRecord
import .Cameras
import .TestScenes

# -- Color constants -- #
BLACK = RGB{Float32}(0.0, 0.0, 0.0)
WHITE = RGB{Float32}(1.0, 1.0, 1.0)

# Ray-Scene intersection:
""" Find the closest intersection point among all objects in the scene
along a ray, constraining the search to values of t between tmin and tmax. """
function closest_intersect(objects::Array{Any, 1}, ray::Ray, tmin, tmax)
    closest_t = tmax
    closest_hit = nothing
    for i in 1:length(objects)
        cur = Scenes.ray_intersect(ray, objects[i])
        if cur !== nothing && cur.t < closest_t && cur.t > tmin
            closest_t = cur.t
            closest_hit = cur
        end
    end 
    return closest_hit
end

""" Approximate gamma correction with gamma = 2 """
function gamma_correct(color)
    r = sqrt(red(color))
    g = sqrt(green(color))
    b = sqrt(blue(color))

    return RGB{Float32}(r,g,b)
end

""" Trace a ray from orig along ray through scene, using Whitted recursive raytracing 
limited to rec_depth recursive calls. """
function traceray(scene::Scene, ray::Ray, tmin, tmax, rec_depth)
    closest_hitrec = closest_intersect(scene.objects, ray, tmin, tmax)

    if closest_hitrec === nothing
        unit_dir = normalize(ray.direction)
        t = 0.5*(unit_dir[2] + 1.0)
        return (1.0-t)*WHITE + t*scene.background
    end
    
    normal = closest_hitrec.normal
    object = closest_hitrec.object
    point = closest_hitrec.intersection
    material = object.material
    
    if rec_depth <= 0
        return BLACK
    end
    color = getColor(scene, material, ray, normal, point, tmax, rec_depth)
    return gamma_correct(color)
end

""" Compute reflection ray direction for reflective surfaces """
function reflectRay(r::Vec3, n::Vec3)
    return 2*dot(n,r)*n - r
end

""" Compute transmission ray direction """
function transmitRay(incident::Vec3, normal::Vec3, index::Float64)
    eta_i = 1.0003 #Refractive index of air 
    eta_t = index
    ref_n = normal
    i_dot_n = dot(incident, normal)

    if (i_dot_n < 0) #Outside surface
        i_dot_n = -i_dot_n
    elseif (i_dot_n > 0) #Inside surface
        ref_n = -normal
        eta_i =  index
        eta_t = 1.0003 #Index for air
    end
    eta = eta_i / eta_t
    discriminant = 1 - (eta*eta)*(1-(i_dot_n*i_dot_n))
    if (discriminant < 0) #Total internal reflection
        return reflectRay(incident, ref_n)
    end
    return (eta * incident) + ((eta * i_dot_n) - sqrt(discriminant)) * normal
end

function getColor(scene::Scene, material::Lambertian, ray::Ray, normal, point, tmax, rec_depth)
    return shade_diffuse(material, ray, normal, point, scene)
end

function getColor(scene::Scene, material::Metallic, ray::Ray, normal, point, tmax, rec_depth)
    # -- Mirror Reflection -- #
    diffuse_color = shade_diffuse(material, ray, normal, point, scene)
    reflectDirection = reflectRay(-ray.direction, normal)
    if (dot(reflectDirection, normal) < 0)
        return BLACK
    end
    reflectionRay = Ray(point, reflectDirection)
    color = traceray(scene, reflectionRay, 1e-8, tmax, rec_depth - 1)
    color = colorMultiply(diffuse_color, color)  
    return color 
end

function getColor(scene::Scene, material::Glossy, ray::Ray, normal, point, tmax, rec_depth)
    reflectDirection = reflectRay(-ray.direction, normal)
    w = normalize(reflectDirection)
    u = normalize(cross(Vec3(0.00424, 1, 0.00764), w))
    v = cross(u, w)
    sp = randomCosineHemi(material.specularExp) #sample hemisphere around reflectDirection
    jittered = sp[1]*u + sp[2]*v + sp[3]*w
    if dot(normal, jittered) < 0
        jittered = (-sp[1])*u + (-sp[2])*v + sp[3]*w
    end
    phongLobe = dot(w, jittered)^material.specularExp
    brdf_sample = material.reflectance*material.specularColor*phongLobe #k_r * c_r * dot(mirror, jittered)^e
    pdf = phongLobe * dot(normal, jittered)
    reflectionRay = Ray(point, jittered)
    color = colorMultiply(brdf_sample, traceray(scene, reflectionRay, 1e-8, tmax, rec_depth - 1))
    return color * (dot(normal, jittered)/pdf)
end

function getColor(scene::Scene, material::Dielectric, ray::Ray, normal, point, tmax, rec_depth)
    theta = acos(dot(-ray.direction, normal)/(norm(ray.direction)*norm(normal)))
    ks = fresnelReflectance(theta, material.index)
    reflectDirection = reflectRay(-ray.direction, normal)
    transmitDirection = transmitRay(ray.direction, normal, material.index)
    reflectionRay = Ray(point, reflectDirection)
    reflectionColor = ks * traceray(scene, reflectionRay, 1e-8, tmax, rec_depth - 1)
    if (transmitDirection == reflectDirection) #TIR case
        return reflectionColor
    else
        transmissionRay = Ray(point, transmitDirection)
        refractionColor = (1-ks)*traceray(scene, transmissionRay, 1e-8, tmax, rec_depth - 1)
    end
    return reflectionColor + refractionColor 
end

""" Compute Fresnel reflectance using Schlick's approximation """ 
function fresnelReflectance(theta::Float64, index::Float64)
    r0 = (1-index) / (1+index)
    r0 *= r0
    return r0 + (1-r0)*(1-cos(theta))^5
end

function glossyReflect(r::Vec3, n::Vec3)
    w_hat = normalize(r)
    #Create non-colinear vector t
    t = [w_hat[1], w_hat[2], w_hat[3]]
    t[argmin(w_hat)] = 1
    t = Vec3(t[1], t[2], t[3])

    u_hat = normalize(cross(t,w_hat))
    v_hat = cross(w_hat, u_hat)
    sampled = randomHemiSphere(n)
    wi = sampled[1]*u_hat + sampled[2]*v_hat + sampled[3]*w_hat
    if dot(n, wi) < 0
        wi = -wi
    end
    return wi 
end

""" Determine the color of interesction point described by hitrec 
Flat shading - just color the pixel the material's diffuse color """
function flat_color(material::Material, hitrec::HitRecord)
    get_diffuse(material, hitrec.uv)
end

""" Normal shading - color-code pixels according to their normals """
function normal_color(normal::Vec3)
    normal_color = normalize(hitrec.normal) / 2 .+ 0.5
    RGB{Float32}(normal_color...)
end

""" Determine if a point is in shadow with a point light """
function is_shadowed(scene, light::PointLight, point::Vec3)
    light_vector = light_direction(light, point)
    shadow_ray = Ray(point, light_vector) #p + tL
    return closest_intersect(scene.objects, shadow_ray, 1e-8, 1) !== nothing
end

""" Determine if a point is in shadow with a directional light """
function is_shadowed(scene, light::DirectionalLight, point::Vec3)
    light_vector = light_direction(light, point)
    shadow_ray = Ray(point, light_vector) #p + tL
    return closest_intersect(scene.objects, shadow_ray, 1e-8, Inf32) !== nothing
end

""" Determine the diffuse color of a physical surface """
function shade_diffuse(material::Material, ray::Ray, normal::Vec3, point::Vec3, scene::Scene)
    color = BLACK 
    for i in 1:length(scene.lights) #for each light in the scene:
        cur_light = scene.lights[i]
        if !is_shadowed(scene, cur_light, point)
            light_val = shade_light(material, ray, normal, point, cur_light) #determine the light's contribution
            color += light_val #add the light's contribution into the color
        end
    end
    return color
end

# Main loop:
function main(scene, camera, height, width, outfile)

    # get the requested scene and camera
    scene = TestScenes.get_scene(scene)
    camera = TestScenes.get_camera(camera, height, width)

    # Create a blank canvas to store the image:
    canvas = zeros(RGB{Float32}, height, width)

    for i in 1:height #for each pixel
        for j in 1:width
            pixel_color = BLACK
            for s in 1:Cameras.SAMPLES_PER_PIXEL #Shoot n rays to each pixel
                subpix = Cameras.convertSubpixel(s)
                view_ray = Cameras.pixel_to_ray(camera, subpix, i, j) #generate viewing ray
                pixel_color += traceray(scene, view_ray, 1, Inf32, 8) #trace the viewing ray to determine pixel color
            end
            pixel_color *= (1/Cameras.SAMPLES_PER_PIXEL) #Average the pixel color across each sample
            canvas[i,j] = pixel_color #set the color of the pixel based on the traced ray
        end
    end

    # clamp canvas to valid range:
    clamp01!(canvas)
    save(outfile, canvas)
end

end # module RayTracer

