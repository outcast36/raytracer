module TestScenes

using ..GfxBase
using ..Scenes
using ..Materials
using ..Lights
using ..MeshGen
using ..Cameras

# -- Color constants -- #
BLACK = RGB{Float32}(0,0,0)
RED = RGB{Float32}(1,0,0)
GREEN = RGB{Float32}(0,1,0)
BLUE = RGB{Float32}(0,0,1)
WHITE = RGB{Float32}(1,1,1)
PURPLE = RGB{Float32}(0.55,0.47,1)
LEMON = RGB{Float32}(1.0, 1.0, 0.3)
PINK = RGB{Float32}(0.7, 0.3, 0.3)
LIME = RGB{Float32}(0.33, 0.97, 0.26)

function camera_1(img_height, img_width)
    CanonicalCamera(img_height, img_width)
end

function camera_2(img_height, img_width)
    eye = Vec3(20, 4, 10)
    view = Vec3(-1, 0, -5) - eye
    up = Vec3(0, 1, 0)
    focal = 8.0
    Cameras.PerspectiveCamera(eye, view, up, focal, img_height, img_width)
end

function camera_3(img_height, img_width)

    Cameras.PerspectiveCamera(
                 Vec3(-1, 0.8, -1.2),  # eye::Vec3
                 Vec3(1, -1, -1), # view::Vec3
                 Vec3(0, 1, 0),   # up::Vec3
                 0.3,     # focal::Real
                 img_height, # canv_height::Int
                 img_width) # canv_width::Int)
end

function camera_4(img_height, img_width)

    Cameras.PerspectiveCamera(
                 Vec3(2.0, 0.4, -0.7),  # eye::Vec3
                 Vec3(-1.0, -0.1, -1), # view::Vec3
                 Vec3(0, 1, 0),   # up::Vec3
                 0.27,     # focal::Real
                 img_height, # canv_height::Int
                 img_width) # canv_width::Int)
end



cameras = [camera_1, camera_2, camera_3, camera_4]

function get_camera(i, img_height, img_width)
    cameras[i](img_height, img_width)
end


function get_scene(i)
    scenes[i]()
end

""" Metallic spheres on the left and right, with a diffuse red sphere in the center """ 
function scene_1()
    bg = RGB{Float32}(0.5, 0.7, 1.0)
    ground_color = RGB{Float32}(0.8, 0.8, 0.0)
    left_metal = RGB{Float32}(0.8, 0.8, 0.8)
    right_metal = RGB{Float32}(0.8, 0.6, 0.2)
    lightColor = RGB{Float32}(0.7, 0.7, 0.7)

    objs = []
    mat_center = Lambertian(RGB{Float32}(0.7, 0.3, 0.3), RED, 10)
    mat_left = Metallic(left_metal, left_metal, 100)
    mat_right = Metallic(right_metal, right_metal, 100)
    mat_ground = Lambertian(ground_color, WHITE, 10)

    push!(objs, Sphere(Vec3(0, 0.0, -3.0), 0.5, mat_center))
    push!(objs, Sphere(Vec3(-1.0, 0.0, -3.5), 0.5, mat_left))
    push!(objs, Sphere(Vec3(1.0, 0.0, -3.5), 0.5, mat_right))
    push!(objs, Sphere(Vec3(0,-101,-1), 100.0, mat_ground))
    lights = [DirectionalLight(lightColor, Vec3(0, 1, 1)),
              PointLight(lightColor, Vec3(1,0,0))]
    Scene(bg, objs, lights)
end

""" Metallic sphere right, with a diffuse blue sphere in the center, and a glass sphere on the left """ 
function scene_2()
    bg = RGB{Float32}(0.5, 0.7, 1.0)
    ground_color = RGB{Float32}(0.8, 0.8, 0.0)
    right_metal = RGB{Float32}(0.8, 0.6, 0.2)
    lightColor = RGB{Float32}(0.6, 0.6, 0.6)

    objs = []
    mat = Lambertian(RGB{Float32}(0.1, 0.2, 0.5), BLUE, 10)
    mat_left = Dielectric(1.52)
    mat_right = Metallic(right_metal, right_metal, 100)
    mat_ground = Lambertian(ground_color, WHITE, 10)

    push!(objs, Sphere(Vec3(0, 0.0, -3.0), 0.5, mat))
    push!(objs, Sphere(Vec3(-1.0, 0.0, -3.5), 0.5, mat_left))
    push!(objs, Sphere(Vec3(1.0, 0.0, -3.5), 0.5, mat_right))
    push!(objs, Sphere(Vec3(0,-101,-1), 100.0, mat_ground))
    lights = [DirectionalLight(lightColor, Vec3(0, 1, 1)),
              PointLight(lightColor, Vec3(1,0,0))]
    Scene(bg, objs, lights)
end

""" Metallic sphere right, with a diffuse blue sphere on the left, and a diamond sphere in the center """ 
function scene_3()
    bg = RGB{Float32}(0.5, 0.7, 1.0)
    ground_color = RGB{Float32}(0.8, 0.8, 0.0)
    right_metal = RGB{Float32}(0.8, 0.6, 0.2)
    lightColor = RGB{Float32}(0.6, 0.6, 0.6)

    objs = []
    mat_center = Dielectric(2.42)
    mat_left = Lambertian(RGB{Float32}(0.1, 0.2, 0.5), BLUE, 10)
    mat_right = Metallic(right_metal, right_metal, 100)
    mat_ground = Lambertian(ground_color, WHITE, 10)

    push!(objs, Sphere(Vec3(0, 0.0, -3.0), 0.5, mat_center))
    push!(objs, Sphere(Vec3(-1.0, 0.0, -3.5), 0.5, mat_left))
    push!(objs, Sphere(Vec3(1.0, 0.0, -3.5), 0.5, mat_right))
    push!(objs, Sphere(Vec3(0,-101,-1), 100.0, mat_ground))
    lights = [DirectionalLight(lightColor, Vec3(0, 1, 1)),
              PointLight(lightColor, Vec3(1,0,0))]
    Scene(bg, objs, lights)
end

""" RGB spheres on the ground with big reflective sphere on the right """
function scene_4()
    bg = BLACK
    lightColor = RGB{Float32}(0.5, 0.5, 0.5)

    r = Lambertian(RED, WHITE, 10)
    g = Lambertian(GREEN, WHITE, 10)
    b = Lambertian(BLUE, WHITE, 10)
    refl = Metallic(WHITE, WHITE, 10)

    objs = [Sphere(Vec3(-1,  -1.1, -3), 0.5, r),
            Sphere(Vec3( -0.5,  -1.0, -4), 0.5, g),
            Sphere(Vec3( 0,  -1.0, -5), 0.5, b),
            Sphere(Vec3( 5,  -1, -4), 4, refl),
            Sphere(Vec3(  0, -101, 0), 100, refl) # floor
           ]

    lights = [DirectionalLight(lightColor, Vec3(-1, 2, 1)),
              PointLight(lightColor, Vec3(-4, -1, 2))]

    Scene(bg, objs, lights)
end


function scene_5()
    bg = RGB{Float32}(0.5, 0.7, 1.0)
    mat_ground = Metallic(WHITE, WHITE, 10)
    exp_100 = Glossy(LIME, LIME, 0.9, 100)
    red_glossy = Glossy(bg, PINK, 0.7, 1000)
    glass = Dielectric(1.53)
    objs = []
    push!(objs, Sphere(Vec3(0,-101,-1), 100.0, mat_ground))
    push!(objs, Sphere(Vec3(1.0, -0.1, -7.5), 1.0, exp_100))
    push!(objs, Sphere(Vec3(0.35, -0.2, -2.9), 0.5, glass))
    push!(objs, Sphere(Vec3(-1.0, -0.1, -6.8), 0.8, red_glossy))
    lights = [DirectionalLight(0.6, Vec3(0, 1, 1)),
    PointLight(0.4, Vec3(1,0,0))]
    Scene(bg, objs, lights)
end

function scene_6()
    bg = RGB{Float32}(0.5, 0.7, 1.0)
    lightColor = RGB{Float32}(0.5, 0.5, 0.5)
    grey = RGB{Float32}(0.8, 0.8, 0.8)
    mat_ground = Metallic(PINK, PINK, 10)
    mat_diamond = Dielectric(2.42)
    exp_10000 = Glossy(bg, grey, 0.65, 10000)
    exp_100 = Glossy(bg, grey, 0.65, 100)

    objs = []
    push!(objs, Sphere(Vec3(0,-101,-1), 100.0, mat_ground))
    push!(objs, Sphere(Vec3(2.0, 0.4, -1.3), 0.15, mat_diamond))
    push!(objs, Sphere(Vec3(-0.4, 0.4, -0.5), 0.5, exp_10000))
    push!(objs, Sphere(Vec3(0.1, 0.45, -1.45), 0.5, exp_100))
    push!(objs, Sphere(Vec3(0.5, 0.5, -2.4), 0.5, exp_100))
    push!(objs, Sphere(Vec3(1.0, 0.55, -3.35), 0.5, exp_10000))
    push!(objs, Sphere(Vec3(1.57, 0.6, -4.4), 0.5, exp_100))

    lights = [DirectionalLight(0.6, Vec3(0, 1, 1)),
    PointLight(0.4, Vec3(1,0,0))]
    Scene(bg, objs, lights)
end

function scene_7()
    bg = RGB{Float32}(0.5, 0.7, 1.0)
    centerColor = RGB{Float32}(0.0, 0.8, 1.0)
    lightColor = RGB{Float32}(0.5, 0.5, 0.5)
    mat_ground = Metallic(WHITE, WHITE, 10)
    mat_center = Lambertian(centerColor, WHITE, 100)

    objs = []
    push!(objs, Sphere(Vec3(0,-101,-1), 100.0, mat_ground))
    push!(objs, Sphere(Vec3(0.0, 0.25, -4.0), 0.6, mat_center))

    lights = [AreaLight(lightColor, Vec3(1, 1, 0), 0.5), DirectionalLight(lightColor, Vec3(0, 1, 1))]
    #lights = [PointLight(lightColor, Vec3(1,1,0))]
    Scene(bg, objs, lights)
end

""" Take the OBJMesh mesh and return an array of Triangles from the mesh
with the given material, after scaling the mesh positions by scale and moving
them by translation """
function mesh_helper(mesh, material, scale=1.0, translation=Vec3(0,0,0))

    for i in 1:length(mesh.positions)
        mesh.positions[i] = mesh.positions[i] * scale + translation
    end

    create_triangles(mesh, material)
end


scenes = [scene_1, scene_2, scene_3, scene_4, scene_5, scene_6, scene_7]

end # module TestScenes
