module Surface

export Sphere, HitRecord, ray_intersect
using LinearAlgebra

using ..GfxBase
using ..MeshGen
using ..Materials

""" Structure to store data about an intersection
of a ray with an object (a "hit")."""
mutable struct HitRecord
    t::Float64
    intersection::Vec3
    normal::Vec3
    object
end

struct Sphere
    center::Vec3
    radius::Float64
    material::Material
end

struct Triangle
    geometry::OBJTriangle
    mesh::OBJMesh
    material::Material
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

""" Intersect a ray with an object.
Returns a HitRecord with info about the intersetion, or nothing if
the ray doesn't intersect with the object. """
function ray_intersect(object, ray::Ray) end

""" Ray-sphere intersection. """
function ray_intersect(object::Sphere, ray::Ray)
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
        surface_normal = normalize(inverse_radius * Vec3(intersection_pt-object.center))
        return HitRecord(t_0, intersection_pt, surface_normal, object)
    end
end

function ray_intersect(object::Triangle, ray::Ray)
    
    return HitRecord(t, intersection, normal, uv, object)
end

end #Module Surface