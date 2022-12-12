module GfxBase

export Vec3, Vec2, Ray, randomUnitSphere, randomHemiSphere

using Images
using StaticArrays
import LinearAlgebra.cross
import LinearAlgebra.normalize
import Random.rand

export RGB

# some type aliases:
const Vec3 = SVector{3, Float64} # 3-vector of floats
const Vec2 = SVector{2, Float64} # 2-vector of floats

struct Ray
    origin::Vec3
    direction::Vec3
end

""" Return a random value in the range [min, max) """
function random_float(low, high)
    return low + (high-low)*rand();
end

""" Rejection method for picking a random point in the unit sphere """
function randomUnitSphere()
    while(true)
        p = Vec3(random_float(-1,1), random_float(-1,1), random_float(-1,1))
        if (norm(p) < 1)
            return p
        end
    end
end

""" Uniformly generate points on the hemisphere around a normal """
function randomHemiSphere(normal::Vec3)
    inSphere = randomUnitSphere()
    if (dot(inSphere, normal) > 0.0) #In the same hemisphere as the normal
        return inSphere;
    else
        return -inSphere;
    end
end

end # module GfxBase
