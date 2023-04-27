module GfxBase

export Vec3, Vec2, Ray, randomUnitSphere, randomHemiSphere, uniformHemiSphere, cosineWeightedHemisphere
export colorMultiply, uniformSampleSphere, randomUnitDisk

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

""" Uniformly sample a point on the hemisphere """
#References: 
#  https://alexanderameye.github.io/notes/sampling-the-hemisphere/
#  https://www.cg.tuwien.ac.at/sites/default/files/course/4411/attachments/04_path_tracing_0.pdf  -- slide 21/115
function uniformHemiSphere()
    sample = Vec2(rand(), rand())
    cosTheta = sample[1]
    sinTheta = sqrt(1-(cosTheta*cosTheta))
    cosPhi = cos(2*pi*sample[2])
    sinPhi = sin(2*pi*sample[2])
    sx = cosPhi*sinTheta
    sy = sinPhi*sinTheta
    sz = cosTheta
    return Vec3(sx, sy, sz)
end

function cosineWeightedHemisphere()
    sample = Vec2(rand(), rand())
    cosTheta = sqrt(sample[1])
    sinTheta = sqrt(1-(cosTheta*cosTheta))
    cosPhi = cos(2*pi*sample[2])
    sinPhi = sin(2*pi*sample[2])
    sx = cosPhi*sinTheta
    sy = sinPhi*sinTheta
    sz = cosTheta
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

function randomUnitDisk()
    while (true)
        p = Vec2(random_float(-1,1), random_float(-1,1))
        if (norm(p) < 1)
            return p
        end
    end
end

end # module GfxBase
