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
SILVER = RGB{Float32}(0.8, 0.8, 0.8)
GOLD = RGB{Float32}(0.8, 0.6, 0.2)
CHERRY = RGB{Float32}(1.0, 0.3, 0.3)

""" Take the OBJMesh mesh and return an array of Triangles from the mesh
with the given material, after scaling the mesh positions by scale and moving
them by translation """
function mesh_helper(mesh, emission, material, scale=1.0, translation=Vec3(0,0,0))

    for i in 1:length(mesh.positions)
        mesh.positions[i] = mesh.positions[i] * scale + translation
    end

    create_triangles(mesh, emission, material)
end

function camera_1()
    Cameras.PerspectiveCamera(
        Vec3(0,0,0), #eye
        Vec3(0,0,-1), #view
        Vec3(0,1,0), #up
        53.1301, #FOV
        1.0, #focal
        16.0, #aperture
        1 #aspect ratio 1:1
    )
end

function camera_2(img_height, img_width)
    eye = Vec3(20, 4, 10)
    view = Vec3(-1, 0, -5) - eye
    up = Vec3(0, 1, 0)
    focal = 8.0
    apt = 16.0
    Cameras.PerspectiveCamera(eye, view, up, focal, apt, img_height, img_width)
end

function camera_3(img_height, img_width)

    Cameras.PerspectiveCamera(
        Vec3(-1, 0.8, -1.2),  # eye::Vec3
        Vec3(1, -1, -1), # view::Vec3
        Vec3(0, 1, 0),   # up::Vec3
        0.3, # focal::Real
        4.0,     
        img_height, # canv_height::Int
        img_width) # canv_width::Int)
end

function camera_4(img_height, img_width)

    Cameras.PerspectiveCamera(
        Vec3(0, 0.0, 0.0),  # eye::Vec3
        Vec3(0.0, 0.0, -1.0), # view::Vec3
        Vec3(0, 1, 0),   # up::Vec3
        4.0,     # focal::Real
        4.0,     #aperture number n
        img_height, # canv_height::Int
        img_width) # canv_width::Int)
end

#function camera_5(img_height, img_width)
#
#    Cameras.OrthographicCamera(
#        Vec3(0, 0.0, 0.0),  # eye::Vec3
#        Vec3(0.0, 0.0, -1.0), # view::Vec3
#        Vec3(0, 1, 0),   # up::Vec3
#        img_height,
#        img_width)
#end

#Parameters sourced from: https://www.graphics.cornell.edu/online/box/data.html
function cornellBoxCam()
    Cameras.PerspectiveCamera(
        Vec3(278, 273, -800), #eye
        Vec3(0,0,1), #view
        Vec3(0,1,0), #up
        39.3076, #VFOV angle in degrees -- roughly a 25mm x 25mm sensor 
        0.035, #focal length (35mm lens)
        4.0, #aperture number
        1, #aspect ratio 1:1
    )
end

cameras = [camera_1, camera_2, camera_3, camera_4, cornellBoxCam]

function get_camera(i)
    cameras[i]()
end


function get_scene(i)
    scenes[i]()
end

function cornellBox()
    bg = BLACK 
    lightColor = RGB{Float32}(0.78, 0.78, 0.78)
    mat_light = Lambertian(lightColor, lightColor, 10)
    objs = []

    #floor
    whiteMat = Lambertian(WHITE, WHITE, 10)
    floor = read_obj("meshes/cb/floor.obj")
    append!(objs, mesh_helper(floor, BLACK, whiteMat))

    #left wall (red)
    redWall = RGB{Float32}(1.0, 0.15, 0.0)
    redMat = Lambertian(redWall, redWall, 10)
    left = read_obj("meshes/cb/left.obj")
    append!(objs, mesh_helper(left, BLACK, redMat))

    #back wall
    back = read_obj("meshes/cb/back.obj")
    append!(objs, mesh_helper(back, BLACK, whiteMat))

    #right wall (green)
    greenWall = RGB{Float32}(0.15, 1.0, 0.0)
    greenMat = Lambertian(greenWall, greenWall, 10)
    right = read_obj("meshes/cb/right.obj")
    append!(objs, mesh_helper(right, BLACK, greenMat))

    #top wall
    top = read_obj("meshes/cb/top.obj")
    append!(objs, mesh_helper(top, BLACK, whiteMat))

    #tall block
    tall = read_obj("meshes/cb/tall_box.obj")
    append!(objs, mesh_helper(tall, BLACK, whiteMat))

    #light source
    lightBox = read_obj("meshes/cb/light.obj")
    append!(objs, mesh_helper(lightBox, WHITE, mat_light))

    #small block
    small = read_obj("meshes/cb/small_box.obj")
    append!(objs, mesh_helper(small, BLACK, whiteMat))

    println("Objects in scene ", length(objs))


    lights = [TriangleLight(f) for f in create_triangles(lightBox, WHITE, mat_light)]
    println("Light objects in scene ", length(lights))
    Scene(bg, objs, lights)
end

function cornellBox2()
    bg = BLACK 
    lightColor = RGB{Float32}(0.78, 0.78, 0.78)
    mat_light = Lambertian(lightColor, lightColor, 10)
    objs = []
    
    reflectMat = Specular(SILVER, SILVER, 10)
    #Spheres
    push!(objs, Sphere(Vec3(160, 110, 380), 110, BLACK, reflectMat))
    push!(objs, Sphere(Vec3(420, 110, 220), 110, BLACK, reflectMat))

    #floor
    whiteMat = Lambertian(WHITE, WHITE, 10)
    floor = read_obj("meshes/cb/floor.obj")
    append!(objs, mesh_helper(floor, BLACK, whiteMat))
    
    #left wall (red)
    redWall = RGB{Float32}(1.0, 0.25, 0.15)
    redMat = Lambertian(redWall, redWall, 10)
    left = read_obj("meshes/cb/left.obj")
    append!(objs, mesh_helper(left, BLACK, redMat))
    
    #back wall
    back = read_obj("meshes/cb/back.obj")
    append!(objs, mesh_helper(back, BLACK, whiteMat))
    
    #right wall (blue)
    blueWall = RGB{Float32}(0.2, 0.2, 1.0)
    blueMat = Lambertian(blueWall, blueWall, 10)
    right = read_obj("meshes/cb/right.obj")
    append!(objs, mesh_helper(right, BLACK, blueMat))
    
    #top wall
    top = read_obj("meshes/cb/top.obj")
    append!(objs, mesh_helper(top, BLACK, whiteMat))
    
    #light source
    lightBox = read_obj("meshes/cb/light.obj")
    append!(objs, mesh_helper(lightBox, WHITE, mat_light))
    
    println("Objects in scene ", length(objs))
    lights = [TriangleLight(f) for f in create_triangles(lightBox, WHITE, mat_light)]
    println("Light objects in scene ", length(lights))
    Scene(bg, objs, lights)
end

""" Red, blue, and green lambertian spheres illuminated by a sphereical light behind the camera """ 
function scene_1()
    bg = BLACK
    #ground_color = RGB{Float32}(0.8, 0.8, 0.0)

    objs = []
    mat_left = Lambertian(RGB{Float32}(0.5, 0.7, 1.0), BLUE, 10) 
    mat_right = Lambertian(LIME, GREEN, 10)
    mat_center = Lambertian(RGB{Float32}(0.7, 0.3, 0.3), RED, 10)
    mat_ground = Lambertian(WHITE, WHITE, 10)
    mat_light = Lambertian(WHITE, WHITE, 10)

    push!(objs, Sphere(Vec3(-1.0, 0.0, -3.5), 0.5, BLACK, mat_left)) #Blue
    push!(objs, Sphere(Vec3(0, 0.0, -3.0), 0.5, BLACK, mat_center)) #Red
    push!(objs, Sphere(Vec3(1.0, 0.0, -3.5), 0.5, BLACK, mat_right)) #Green
    push!(objs, Sphere(Vec3(0,-101,-1), 100.0, BLACK, mat_ground))
    push!(objs, Sphere(Vec3(0,3.0,-5.0), 1.0, WHITE, mat_light)) #Light source

    lights = [SphericalLight(objs[5])] 
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

    push!(objs, Sphere(Vec3(0, 0.0, -3.0), 0.5, mat, BLACK))
    push!(objs, Sphere(Vec3(-1.0, 0.0, -3.5), 0.5, mat_left, BLACK))
    push!(objs, Sphere(Vec3(1.0, 0.0, -3.5), 0.5, mat_right, BLACK))
    push!(objs, Sphere(Vec3(0,-101,-1), 100.0, mat_ground, BLACK))
    push!(objs, Sphere(Vec3(1.5, 1.5, 0.0), 2.0, nothing, BLACK))

    lights = [AreaLight(lightColor, )]
    #lights = [DirectionalLight(lightColor, Vec3(0, 1, 1)), PointLight(lightColor, Vec3(1,0,0))]
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

    lights = [AreaLight(lightColor, Vec3(1.5, 1.5, 0.0), 2.0)]
    #lights = [DirectionalLight(lightColor, Vec3(0, 1, 1)), PointLight(lightColor, Vec3(1,0,0))]
    Scene(bg, objs, lights)
end

""" RGB spheres on the ground with big reflective sphere on the right """
function scene_4()
    bg = BLACK
    lb = RGB{Float32}(0.5, 0.7, 1.0)

    r = Lambertian(RED, WHITE, 10)
    g = Lambertian(GREEN, WHITE, 10)
    b = Lambertian(BLUE, WHITE, 10)
    refl = Lambertian(WHITE, WHITE, 10)#Metallic(WHITE, WHITE, 10)

    objs = [Sphere(Vec3(-1,  -1.1, -3), 0.5, BLACK, r),
            Sphere(Vec3( -0.5,  -1.0, -4), 0.5, BLACK, g),
            Sphere(Vec3( 0,  -1.0, -5), 0.5, BLACK, b),
            Sphere(Vec3( 5,  -1, -4), 4, BLACK, refl),
            Sphere(Vec3(  0, -101, 0), 100, BLACK, refl), 
            Sphere(Vec3(0,0,0), 0.8, WHITE, refl), #light
           ]

    lights = [AreaLight(objs[6])]
    Scene(bg, objs, lights)
end


function scene_5()
    bg = RGB{Float32}(0.5, 0.7, 1.0)
    lightColor = RGB{Float32}(0.5, 0.5, 0.5)
    mat_ground = Metallic(WHITE, WHITE, 10)
    exp_100 = Glossy(LIME, LIME, 0.9, 100)
    red_glossy = Glossy(bg, PINK, 0.7, 1000)
    glass = Dielectric(1.53)
    objs = []
    push!(objs, Sphere(Vec3(0,-101,-1), 100.0, mat_ground))
    push!(objs, Sphere(Vec3(1.0, -0.1, -7.5), 1.0, exp_100))
    push!(objs, Sphere(Vec3(0.35, -0.2, -2.9), 0.5, glass))
    push!(objs, Sphere(Vec3(-1.0, -0.1, -6.8), 0.8, red_glossy))
    
    lights = [AreaLight(lightColor, Vec3(1.5, 1.5, 0.0), 2.0)]
    #lights = [DirectionalLight(0.6, Vec3(0, 1, 1)), PointLight(0.4, Vec3(1,0,0))]
    Scene(bg, objs, lights)
end

function scene_6()
    bg = RGB{Float32}(0.5, 0.7, 1.0)
    centerColor = RGB{Float32}(0.0, 0.8, 1.0)
    lightColor = RGB{Float32}(0.6, 0.6, 0.6)
    bronze = RGB{Float32}(0.72, 0.45, 0.2)
    mat_ground = Metallic(WHITE, WHITE, 10)
    mat_blue = Lambertian(centerColor, WHITE, 100)
    mat_diamond = Dielectric(2.42)
    mat_bronze = Metallic(bronze, bronze, 1000)

    objs = []
    push!(objs, Sphere(Vec3(0,-101,-1), 100.0, mat_ground))
    push!(objs, Sphere(Vec3(1.0, -0.1, -3.2), 0.5, mat_bronze))
    push!(objs, Sphere(Vec3(0.2, 0.1, -4.5), 0.7, mat_diamond))
    push!(objs, Sphere(Vec3(-1.1, 0.3, -6.5), 1.1, mat_blue))

    lights = [AreaLight(lightColor, Vec3(1.5, 1.5, 0.0), 2.0)] #
    #lights = [PointLight(lightColor, Vec3(-0.025, 2.033, -1.178)), DirectionalLight(lightColor, Vec3(0, 1, 1))]
    Scene(bg, objs, lights)
end

function scene_7()
    bg = BLACK
    lightColor = RGB{Float32}(0.6, 0.6, 0.6)
    objs = []

    # add a bunny:
    bunnyColor = RGB{Float32}(0.6, 0.5, 0.5)
    bunny_mat = Lambertian(bunnyColor, bunnyColor, 10)
    bunny = read_obj("meshes/bunny2.obj")
    append!(objs, mesh_helper(bunny, bunny_mat, 1.0, Vec3(0.2, 0, -5)))

    # add a cube
    cube_mat = Metallic(WHITE, WHITE, 10)
    append!(objs, mesh_helper(cube_mesh(), cube_mat, 10.0, Vec3(-11.2, 0, 0)))

    #lights = [AreaLight(lightColor, Vec3(1.5, 1.5, 0.0), 2.0)] #

    lights = [ PointLight(0.5, Vec3(1,2,-5)),
               DirectionalLight(0.3, Vec3(0,0,1)),
               DirectionalLight(0.3, Vec3(0,1,1)),
               DirectionalLight(0.3, Vec3(1,1,1)),
               DirectionalLight(0.3, Vec3(0,1,0)) ]


    Scene(bg, objs, lights)

end


scenes = [scene_1, scene_2, scene_3, scene_4, scene_5, scene_6, scene_7, cornellBox, cornellBox2]

end # module TestScenes
