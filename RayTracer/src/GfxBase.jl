module GfxBase

export Vec3, Vec2, Ray

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

end # module GfxBase
