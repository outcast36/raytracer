module Scenes

export HitRecord, Sphere, Scene, TriangleMesh, ray_intersect, create_triangles
using LinearAlgebra

using ..GfxBase
using ..MeshGen
using ..Materials

""" Scene datatype to store hittable objects, and light sources """
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

struct Sphere
    center::Vec3
    radius::Float64
    emission::RGB{Float32} #If surface is a light emitting surface (a luminaire) this will be nonzero
    material
end

""" Data type: stores the OBJTriangle, a reference to its Mesh
object, and the material it should be rendered with. """
struct Triangle
    geometry::OBJTriangle
    mesh::OBJMesh
    emission::RGB{Float32} #If surface is a light emitting surface (a luminaire) this will be nonzero
    material
end

""" Return an Array of all Triangles belonging to the given mesh, assigning
each one the given material. """
function create_triangles(mesh::OBJMesh, material)
    [Triangle(f, mesh, RGB{Float32}(0,0,0), material) for f in mesh.triangles]
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

""" Intersect a ray with an object.
Returns a HitRecord with info about the intersetion, or nothing if
the ray doesn't intersect with the object. """
function ray_intersect(ray::Ray, object) end

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

        ##Textures
        #d = normalize(intersection_pt- object.center) #D is vector from intersection point to sphere center
        #u = 0.5 + (atan(d[1], d[3])/(2*pi))
        #v = 0.5 + (-asin(d[2])/pi)
        #uv = Vec2(v,u)
        return HitRecord(t_0, intersection_pt, surface_normal, nothing, object)
    end
end

""" Ray-triangle intersection via barycentric coordinates """
function ray_intersect(ray::Ray, object::Triangle)
    #Get vertices
    a = get_vertex(object,1)
    b = get_vertex(object,2)
    c = get_vertex(object,3)

    col_u = a-b
    col_v = a-c
    col_b = a - ray.origin
    #Solve linear system Ax = b via Cramer's rule. Where A is a 3x3 matrix with columns (a-b), (a-c), and d. 
    #x is a column vector: [beta, gamma, t]^T, and b is the column vector (a-p)
    ei_minus_hf = (col_v[2]*ray.direction[3]) - (ray.direction[2]*col_v[3])
    gf_minus_di = (ray.direction[1]*col_v[3]) - (col_v[1]*ray.direction[3])
    dh_minus_eg = (col_v[1]*ray.direction[2]) - (col_v[2]*ray.direction[1])

    ak_minus_jb = (col_u[1]*col_b[2]) - (col_b[1]*col_u[2])
    jc_minus_al = (col_b[1]*col_u[3]) - (col_u[1]*col_b[3])
    bl_minus_kc = (col_u[2]*col_b[3]) - (col_b[2]*col_u[3])

    M = (col_u[1]*ei_minus_hf) + (col_u[2]*gf_minus_di) + (col_u[3]*dh_minus_eg) #Determinant of matrix A 
    beta = ((col_b[1]*ei_minus_hf)+(col_b[2]*gf_minus_di)+(col_b[3]*dh_minus_eg))/M
    gamma = ((ray.direction[3]*ak_minus_jb)+(ray.direction[2]*jc_minus_al)+(ray.direction[1]*bl_minus_kc))/M
    t = -((col_v[3]*ak_minus_jb)+(col_v[2]*jc_minus_al)+(col_v[1]*bl_minus_kc))/M

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
    return HitRecord(t, intersection, normal, nothing, object)
end

end # module Scenes
