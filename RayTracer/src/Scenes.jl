module Scenes

export HitRecord, Sphere, Scene, TriangleMesh, ray_intersect, create_triangles
#export has_uvs, has_normals, get_vertex, get_uv, get_normal

using LinearAlgebra

#push!(LOAD_PATH, pwd())
using ..GfxBase
using ..MeshGen
using ..Materials



#####################################
###### Generic Scene Data Type ######
#####################################
struct Scene
    background::RGB{Float32}
    objects::Array{Any,1}
    lights::Array{Any,1}
end

""" Structure to store data about an intersection
of a ray with an object (a "hit")."""
mutable struct HitRecord
    t::Float64
    intersection::Vec3
    normal::Vec3
    uv::Union{Vec2,Nothing}
    object
end

# Abstract ray-object intersection function:
# Needs to be implemented for each type of object to be rendered
""" Intersect a ray with an object.
Returns a HitRecord with info about the intersetion, or nothing if
the ray doesn't intersect with the object. """
function ray_intersect(ray::Ray, object) end


##################
##### Sphere #####
##################

# Data type:
struct Sphere
    center::Vec3
    radius::Float64
    material::Material
end

""" Ray-sphere intersection. """
function ray_intersect(ray::Ray, object::Sphere)
    squared_radius = object.radius * object.radius
    a = norm(ray.direction)*norm(ray.direction)
    half_b = dot(ray.direction, (ray.origin - object.center))
    c = -squared_radius + norm(ray.origin - object.center)*norm(ray.origin - object.center)
    discriminant = (half_b*half_b) - c*a
    inverse_radius = 1/object.radius
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
        surface_normal = inverse_radius * Vec3(intersection_pt-object.center)

        #Textures
        d = normalize(intersection_pt- object.center) #D is vector from intersection point to sphere center
        u = 0.5 + (atan(d[1], d[3])/(2*pi))
        v = 0.5 + (-asin(d[2])/pi)
        uv = Vec2(v,u)
        return HitRecord(t_0, intersection_pt, surface_normal, uv, object)
    end
end


###########################
###### Triangle Mesh ######
###########################

""" Data type: stores the OBJTriangle, a reference to its Mesh
object, and the material it should be rendered with. """
struct Triangle
    geometry::OBJTriangle
    mesh::OBJMesh
    material
end

""" Return an Array of all Triangles belonging to the given mesh, assigning
each one the given material. """
function create_triangles(mesh::OBJMesh, material)
    [Triangle(f, mesh, material) for f in mesh.triangles]
end

""" Some helper functions that make for easier access to triangle data: """
function get_vertex(tri::Triangle, i)
    tri.mesh.positions[tri.geometry.positions[i]]
end

function has_uvs(tri::Triangle)
    length(tri.geometry.uvs) == 3
end

function get_uv(tri::Triangle, i)
    tri.mesh.uvs[tri.geometry.uvs[i]]
end

function has_normals(tri::Triangle)
    length(tri.geometry.normals) == 3
end

function get_normal(tri::Triangle, i)
    tri.mesh.normals[tri.geometry.normals[i]]
end


function ray_intersect(ray::Ray, object::Triangle)
    #Get vertices
    a = get_vertex(object,1)
    b = get_vertex(object,2)
    c = get_vertex(object,3)

    vec_b = a - ray.origin
    #vec_b = reshape(vec_b, length(vec_b), 1) #Set up the right-hand side 3-vector: b = a-p
    A = hcat((a-c), (a-b), ray.direction) #Set up the matrix A with columns a - b, a - c, and d -- something to do with the ordering of vertices?
    x = inv(A) * vec_b #Why does this work but not the line below
    #x = A / vec_b #Solve the system Ax = b
    #x vector: [beta, gamma, t]^T
    beta = x[1]
    gamma = x[2]
    t = x[3]

    point_in_triangle = (beta > 0 && gamma > 0) && (beta + gamma < 1)
    if !point_in_triangle || t <= 0 
        return nothing
    end
    alpha = (1 - (beta + gamma))
    if has_normals(object)
        normal = alpha * get_normal(object,1) + beta * get_normal(object,2) + gamma * get_normal(object,3)
    else
       normal = normalize(cross((b-a), (c-a)))
    end
    intersection = ray.origin + t * ray.direction

    #Textures
    if has_uvs(object)
        uv =  alpha* get_uv(object,1) + beta * get_uv(object,3) + gamma * get_uv(object,2)
        uv = Vec2(1-uv[2],uv[1])
    else
        uv = nothing
    end
    return HitRecord(t, intersection, normal, uv, object)
end

end # module Scenes
