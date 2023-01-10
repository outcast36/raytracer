module GfxBase

export Vec3, Vec2, Ray, randomUnitSphere, randomHemiSphere, randomCosineHemi, colorMultiply, uniformSampleSphere

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
    return low + (high-low)*rand()
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

""" Sample a point on the hemisphere with cosine density """
function randomCosineHemi(exp)
    sample = Vec2(rand(), rand())
    cos_phi = cos(2*pi*sample[1])
    sin_phi = sin(2*pi*sample[1])
    cos_theta = (1-sample[2])^(1/(1+exp))
    sin_theta = sqrt(1 - (cos_theta*cos_theta))
    sx = sin_theta*cos_phi
    sy = sin_theta*sin_phi
    sz = cos_theta
    return Vec3(sx, sy, sz)
end

""" Multiply two RGB colors channelwise """
function colorMultiply(a, b)
    return RGB{Float32}(a.r*b.r, a.g*b.g, a.b*b.b)
end

""" Uniformly sample over sphere's surface area """
function uniformSampleSphere(sp::Vec2)
    z = 1 - 2*sp[1]
    r = sqrt(max(0, 1 - z^2))
    phi = 2 * pi * sp[2]
    return Vec3(r*cos(phi), r*sin(phi), z)
end

end # module GfxBase
