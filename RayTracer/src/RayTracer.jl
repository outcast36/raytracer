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
import Random.rand

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
        if ((cur !== nothing) && (cur.t < closest_t) && (cur.t > tmin))
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

""" Core of Monte Carlo Path Tracing -- Integrate rendering equation via MC methods """
function traceray(scene::Scene, ray::Ray, throughput, tmin, tmax, depth)
    probrr = max(red(throughput), green(throughput), blue(throughput)) #probability that a path continues bouncing
    #Trace ray to find point of intersection with the nearest surface
    closestHit = closest_intersect(scene.objects, ray, tmin, tmax)
    if (closestHit === nothing)
        return scene.background #Black background
    end
    normal = closestHit.normal
    object = closestHit.object
    point = closestHit.intersection
    material = object.material
    color = object.emission #Constant color in place of emission function (Le)

    #Russian roulette -- Randomly decide to compute emitted or reflected light
    #https://graphics.stanford.edu/courses/cs348b-01/course29.hanrahan.pdf
    roulette = rand()
    termCond = ((roulette > probrr) && (depth >=4)) # (1 - p_rr) probability to terminate path at current step 
    if (termCond)
        return colorMultiply(color, throughput) #If emitted: return weight * Le (Page 9, Step 3A)
    end
    brdfVals = sample(material, -ray.direction, normal)
    brdfFactor = brdfVals.mult
    omega = brdfVals.omega
    throughput = colorMultiply(throughput, brdfFactor) * (1/probrr) #If reflected: weight *= reflectance (Page 9, Step 3B)
    recurRay = Ray(point, omega)
    recurColor = traceray(scene, recurRay, throughput, tmin, tmax, depth+1)
    color += colorMultiply(brdfFactor, recurColor) * (1/probrr) #fr * cos(omega) * 1/P(omega) 
    return color 
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

# Main loop:
function main(scene, camera, width, spp, outfile)

    # get the requested scene and camera
    scene = TestScenes.get_scene(scene)
    camera = TestScenes.get_camera(camera)

    # Create a blank canvas to store the image:
    height = div(width, camera.aspectRatio)
    canvas = zeros(RGB{Float32}, height, width)
    n = isqrt(spp)

    for i in 1:height #for each pixel
        for j in 1:width
            pixel_color = BLACK
            # n x n samples per pixel
            for si in 1:n 
                for sj in 1:n #Loop over n x n subpixel grid
                #Choose a ray given (x, y, u, v, t)
                # x, y: world coordinates of pixel sample for antialiasing
                # u, v: camera aperture
                # t: time of ray for motion blur
                    view_ray = Cameras.pixel_to_ray(camera, height, width, n, si, sj, i, j)
                    pixel_color += traceray(scene, view_ray, WHITE, 1e-8, Inf32, 1)
                end
            end
            pixel_color *= (1/spp)
            canvas[i,j] = gamma_correct(pixel_color) #set the color of the pixel based on the traced ray
        end
    end

    # clamp canvas to valid range:
    clamp01nan!(canvas)
    save(outfile, canvas)
end

end # module RayTracer