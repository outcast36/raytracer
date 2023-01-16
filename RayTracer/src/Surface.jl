module Surface

export Sphere, HitRecord, ray_intersect
using LinearAlgebra

""" Structure to store data about an intersection
of a ray with an object (a "hit")."""
mutable struct HitRecord
    t::Float64
    intersection::Vec3
    normal::Vec3
    object
end

""" Sphere object """
struct Sphere
    center::Vec3
    radius::Float64
    material::Material
end

# Abstract ray-object intersection function:
# Needs to be implemented for each type of object to be rendered
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

end