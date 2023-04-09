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
export pixel_to_ray

""" The most basic perspective camera. The eye is at (0, 0, 0),
the up vector is (0, 1, 0), and the view direction is (0, 0, -1).
The viewport has width 1 and whatever height makes pixels square.
The focal length (distance from eye to viewport) is 1. """


""" A general perspective camera."""
mutable struct PerspectiveCamera

    # orthonormal basis specifies position and orientation:
    eye::Vec3  # position of the eye
    u_axis::Vec3 # points right in image space
    v_axis::Vec3 # points up in image space
    w_axis::Vec3 # opposite the view direction

    vfov::Real # vertical field of view in radians

    # the viewport parallel to the uv plane and centered at (0, 0, -d): -- A2 extension support viewports not centered on this point
    # also support oblique viewing with an additional image plane normal parameter
    # d is distance from camera eye to viewport (image plane) -- No explicit parameter treat image plane as the focal plane

    #Lens diameter is F/n: where F is focal length and n is aperture number (Distributed Ray Tracing, Cook, 1984)
    focalLength::Real #The distance from the lens to the focal plane (f).
    apertureNumber::Real #Ratio between focal length and diameter of the entrance pupil -- same size as lens diameter above

    # image dimensions:
    aspectRatio::Real #Aspect ratio of the viewport and final image
end

#""" An orthographic camera where rays all travel parallel to the view ray, but originate from the center of each (sub)pixel."""
#mutable struct OrthographicCamera
#
#    # orthonormal basis specifies position and orientation:
#    eye::Vec3  # position of the eye
#    u_axis::Vec3 # points right in image space
#    v_axis::Vec3 # points up in image space
#    w_axis::Vec3 # opposite the view direction
#
#    # image dimensions:
#    canv_height::Int # image height in pixels
#    canv_width::Int # image width in pixels
#end

""" Constructors for cameras """
function PerspectiveCamera(eye::Vec3, view::Vec3, up::Vec3, vfov::Real, focal::Real, apertureNumber::Real, aspectRatio::Real)
    w_axis = normalize(-view)
    u_axis = normalize(cross(up, w_axis))
    v_axis = cross(w_axis, u_axis)
    vfov = deg2rad(vfov)
    return PerspectiveCamera(eye, u_axis, v_axis, w_axis, vfov, focal, apertureNumber, aspectRatio)
end

#function OrthographicCamera(eye::Vec3, view::Vec3, up::Vec3, canv_height::Int, canv_width::Int)
#    w_axis = normalize(-view)
#    u_axis = normalize(cross(up, w_axis))
#    v_axis = normalize(cross(w_axis, u_axis))
#    return OrthographicCamera(eye, u_axis, v_axis, w_axis, canv_height, canv_width)
#end

function pixel_to_ray(camera::PerspectiveCamera, imgHeight::Int, imgWidth::Int, n, si, sj, i, j)
    # viewport height = 2 * imagePlaneDist * tan(VFOV/2) -- in this case imagePaneDist == focalLength
    vp_height = 2 * camera.focalLength * tan(camera.vfov/2)
    vp_width = vp_height * camera.aspectRatio

    u = ((j - (sj - rand())/n)/imgWidth - 0.5) * vp_width 
    v = -((i - (si - rand())/n)/imgHeight - 0.5) * vp_height 

    lensPoint = camera.focalLength/(2*camera.apertureNumber) * randomUnitDisk() #lens radius is F/n * 0.5
    px = u#*camera.focalLength
    py = v#*camera.focalLength
    ray_dir = (px - lensPoint[1])*camera.u_axis + (py - lensPoint[2])*camera.v_axis - camera.focalLength*camera.w_axis
    #ray_dir = u*camera.u_axis + v*camera.v_axis - camera.focalLength*camera.w_axis
    ray_origin = camera.eye + lensPoint[1]*camera.u_axis + lensPoint[2]*camera.v_axis #camera.eye + offset
    return Ray(ray_origin, ray_dir)    
end

#function pixel_to_ray(camera::OrthographicCamera, subPixel::Vec2, i, j)
#    # viewport height = 1
#    vp_height = camera.canv_height / camera.canv_width
#
#    # convert from i,j 2D array indices to (u, v) coordinates
#    # with (0,0) in the center of the image
#    # half pixel shift, scale down to size 1, shift origin to center
#    u = ((j - subPixel[2]/N_SUBPIXELS) / (camera.canv_width) - 0.5)
#
#    # same as u, except i increases downward so flip with a negative sign
#    # and multiply by vp_height to account for aspect ratio
#    v = -((i - subPixel[1]/N_SUBPIXELS) / (camera.canv_height) - 0.5) * vp_height
#
#    ray_origin = camera.eye + u*camera.u_axis + v*camera.v_axis #camera.eye + offset
#    return Ray(ray_origin, -camera.w_axis)
#end


end # module Cameras