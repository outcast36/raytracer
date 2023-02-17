""" The Cameras module is responsible for ray generation. It defines
types for each kind of camera and implements a pixel_to_ray function
that takes a camera and a set of pixel indices and returns a ray
along which to search for objects that appear at that pixel.
"""
module Cameras

import LinearAlgebra.normalize
import LinearAlgebra.norm
import LinearAlgebra.cross

#push!(LOAD_PATH, pwd())
using ..GfxBase

export CanonicalCamera
export pixel_to_ray, convertSubpixel
export SAMPLES_PER_PIXEL

#Divide pixel into n*n grid of subpixels
SAMPLES_PER_PIXEL = 256
N_SUBPIXELS = isqrt(SAMPLES_PER_PIXEL)

""" The most basic perspective camera. The eye is at (0, 0, 0),
the up vector is (0, 1, 0), and the view direction is (0, 0, -1).
The viewport has width 1 and whatever height makes pixels square.
The focal length (distance from eye to viewport) is 1. """
mutable struct CanonicalCamera
    canv_height::Int
    canv_width::Int
end

""" Convert a 1d subpixel index into a 2d subpixel index """
function convertSubpixel(subpix)
    sub_i = subpix%N_SUBPIXELS == 0 ? div(subpix,N_SUBPIXELS) : div(subpix,N_SUBPIXELS)+1
    sub_j = subpix%N_SUBPIXELS == 0 ? N_SUBPIXELS : subpix%N_SUBPIXELS 
    return Vec2(N_SUBPIXELS-sub_i, N_SUBPIXELS-sub_j)
end

""" Given a camera and pixel coordinates 1-indexed, i=row,j=column],
return the ray along which to search for objects that project to
that pixel. """
function pixel_to_ray(camera::CanonicalCamera, subPixel::Vec2, i, j)
    # viewport height = 1
    vp_height = camera.canv_height / camera.canv_width

    # convert from i,j 2D array indices to (u, v) coordinates
    # with (0,0) in the center of the image
    # half pixel shift, scale down to size 1, shift origin to center
    u =  ((j - subPixel[2]/N_SUBPIXELS) / (camera.canv_width) - 0.5)

    # same as u, except i increases downward so flip with a negative sign
    # and multiply by vp_height to account for aspect ratio
    v = -((i - subPixel[1]/N_SUBPIXELS) / (camera.canv_height) - 0.5) * vp_height

    # focal length is 1, pointing down the -z axis
    w = -1.0

    Ray(Vec3(0, 0, 0), Vec3(u, v, w))
end

""" A general perspective camera."""
mutable struct PerspectiveCamera

    # orthonormal basis specifies position and orientation:
    eye::Vec3  # position of the eye
    u_axis::Vec3 # points right in image space
    v_axis::Vec3 # points up in image space
    w_axis::Vec3 # opposite the view direction

    # the viewport is width 1, parallel to the uv plane and
    # centered at (0, 0, -focal):
    focalLength::Real #The distance from the lens to the focal plane.
    aperture_n::Real #Aperture number n

    # image dimensions:
    canv_height::Int # image height in pixels
    canv_width::Int # image width in pixels
end

""" An orthographic camera where rays all travel parallel to the view ray, but originate from the center of each (sub)pixel."""
mutable struct OrthographicCamera

    # orthonormal basis specifies position and orientation:
    eye::Vec3  # position of the eye
    u_axis::Vec3 # points right in image space
    v_axis::Vec3 # points up in image space
    w_axis::Vec3 # opposite the view direction

    # image dimensions:
    canv_height::Int # image height in pixels
    canv_width::Int # image width in pixels
end

""" Constructors for cameras """
function PerspectiveCamera(eye::Vec3, view::Vec3, up::Vec3, focal::Real, aperture_n::Real, canv_height::Int, canv_width::Int)
    w_axis = normalize(-view)
    u_axis = normalize(cross(up, w_axis))
    v_axis = normalize(cross(w_axis, u_axis))
    return PerspectiveCamera(eye, u_axis, v_axis, w_axis, focal, aperture_n, canv_height, canv_width)
end

function OrthographicCamera(eye::Vec3, view::Vec3, up::Vec3, canv_height::Int, canv_width::Int)
    w_axis = normalize(-view)
    u_axis = normalize(cross(up, w_axis))
    v_axis = normalize(cross(w_axis, u_axis))
    return OrthographicCamera(eye, u_axis, v_axis, w_axis, canv_height, canv_width)
end


function pixel_to_ray(camera::PerspectiveCamera, subPixel::Vec2, i, j)
    # viewport height = 1
    vp_height = camera.canv_height / camera.canv_width

    # convert from i,j 2D array indices to (u, v) coordinates
    # with (0,0) in the center of the image
    # half pixel shift, scale down to size 1, shift origin to center
    u = ((j - subPixel[2]/N_SUBPIXELS) / (camera.canv_width) - 0.5)

    # same as u, except i increases downward so flip with a negative sign
    # and multiply by vp_height to account for aspect ratio
    v = -((i - subPixel[1]/N_SUBPIXELS) / (camera.canv_height) - 0.5) * vp_height

    lensPoint = camera.focalLength/(2*camera.aperture_n) * randomUnitDisk() #lens diameter is F/n where F is focal length of the lens and n is the aperture number (Cook 1984).
    px = u*(camera.focalLength)
    py = v*(camera.focalLength) #/d?
    ray_dir = (px - lensPoint[1])*camera.u_axis + (py - lensPoint[2])*camera.v_axis - camera.focalLength*camera.w_axis
    ray_dir = normalize(ray_dir)
    ray_origin = camera.eye + lensPoint[1]*camera.u_axis + lensPoint[2]*camera.v_axis #camera.eye + offset
    #ray_dir = u*camera.u_axis + v*camera.v_axis - camera.focalLength*camera.w_axis
    return Ray(ray_origin, ray_dir)
end

function pixel_to_ray(camera::OrthographicCamera, subPixel::Vec2, i, j)
    # viewport height = 1
    vp_height = camera.canv_height / camera.canv_width

    # convert from i,j 2D array indices to (u, v) coordinates
    # with (0,0) in the center of the image
    # half pixel shift, scale down to size 1, shift origin to center
    u = ((j - subPixel[2]/N_SUBPIXELS) / (camera.canv_width) - 0.5)

    # same as u, except i increases downward so flip with a negative sign
    # and multiply by vp_height to account for aspect ratio
    v = -((i - subPixel[1]/N_SUBPIXELS) / (camera.canv_height) - 0.5) * vp_height

    ray_origin = normalize(camera.eye + lensPoint[1]*camera.u_axis + lensPoint[2]*camera.v_axis) #camera.eye + offset
    return Ray(ray_origin, -camera.w_axis)
end


end # module Cameras