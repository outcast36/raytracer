module MeshGen

export read_obj, write_obj
export gen_mesh, est_normals
export cube_mesh, cylinder_mesh, sphere_mesh, estimate_normals
export OBJTriangle, OBJMesh

using FileIO
using LinearAlgebra

push!(LOAD_PATH, pwd())

include("GfxBase.jl")
using .GfxBase


""" OBJTriangle
A struct that represents a single triangle in a mesh. """
mutable struct OBJTriangle
    positions::Array{Int, 1} # vertex position indices
    uvs::Array{Int, 1} # vertex texture coordinate indices
    normals::Array{Int, 1} # normal vector indices
end

""" OBJMesh
A struct that represents an indexed triangle mesh for reading from
or writing to OBJ format. """
mutable struct OBJMesh
    positions::Array{Vec3, 1} # all vertex positions
    uvs::Array{Vec2, 1} # all texture coordinates
    normals::Array{Vec3, 1} # all vertex normals
    triangles::Array{OBJTriangle, 1} # the OBJTriangles belonging to the mesh
end

""" read_obj(obj_filename)
Read a mesh in OBJ format from file obj_filename."""
function read_obj(obj_filename)
    m = OBJMesh([], [], [], []) # create a mesh
    open(obj_filename) do f
        for (line_number, line) in enumerate(eachline(f))
            if line == "" || line[1] == "#"
                continue # skip comments
            end
            # Read the line and add its contents to the correct field of m:
            tokens = split(strip(line))
            if tokens[1] == "v" # vertex
                push!(m.positions, Vec3([parse(Float64, x) for x in tokens[2:end]]...))
            elseif tokens[1] == "vt" # vertex texture
                push!(m.uvs, Vec2([parse(Float64, x) for x in tokens[2:end]]...))
            elseif tokens[1] == "vn" # vertex normal
                push!(m.normals, Vec3([parse(Float64, x) for x in tokens[2:end]]...))
            elseif tokens[1] == "f"
                # create a OBJTriangle face:
                points = []
                uvs = []
                normals = []
                # handle faces with no texture and/or normals
                for corner in tokens[2:end]
                    indices = split(corner, '/')
                    if length(indices) == 3 # all 3 present, third is a normal
                        push!(normals, parse(Int, indices[3]))
                    end
                    if length(indices) >= 2 && indices[2] != ""
                        # if there are 2 or more and the second isn't blank, it's a texture
                        push!(uvs, parse(Int, indices[2]))
                    end
                    if length(indices) >= 1 # first value is the position
                        push!(points, parse(Int, indices[1]))
                    else # unless it has none, in which case it's not valid
                        error("in line $line_number: face vertex $corner could not be parsed")
                    end
                end
                # create the triangle and add it to the triangles array
                push!(m.triangles, OBJTriangle(points, uvs, normals))
            end
        end
    end
    return m
end

""" write_obj(obj_filename)
Write the given mesh in OBJ format to file obj_filename."""
function write_obj(obj_filename, mesh::OBJMesh)
    open(obj_filename, "w") do f
        # write all positions:
        for v in mesh.positions
            write(f, "v $(v[1]) $(v[2]) $(v[3])\n")
        end

        # write all texture coords:
        for v in mesh.uvs
            write(f, "vt $(v[1]) $(v[2])\n")
        end
        # write all normals:
        for v in mesh.normals
            write(f, "vn $(v[1]) $(v[2]) $(v[3])\n")
        end

        # write all triangles:
        for tri in mesh.triangles
            write(f, "f $(tri_vertex_str(tri))\n")
        end

    end

end

""" tri_vertex_str(triangle)
Return a string with the indices of applicable positions, texture coordinates,
and normals for a given triangle according to the OBJ specification.
In particular, if p, u, and n are position, vertex and normal, each corner
of the triangle is represented as one of the following:
    p       (position only)
    p/u     (position and texture)
    p//n    (position and normal)
    p/u/n   (position, texture, and normal) """
function tri_vertex_str(triangle::OBJTriangle)
    # determine whether textures and normals are present:
    write_uv = length(triangle.uvs) == length(triangle.positions)
    write_normals = length(triangle.normals) == length(triangle.positions)
    corners = []
    for i = 1:3
        output = "$(triangle.positions[i])"
        if write_uv && !write_normals
            output = output * "/$(triangle.uvs[i])" # * does concatenation(!)
        elseif !write_uv && write_normals
            output = output * "//$(triangle.normals[i])"
        elseif write_uv && write_normals
            output = output * "/$(triangle.uvs[i])/$(triangle.normals[i])"
        end
        push!(corners, output)
    end
    join(corners, " ")
end


""" gen_mesh(outfile, geom, divisionsU, divisionsV)
Generate a mesh and save the result in a file with name outfile.
geom may be "cube", "cylinder", or "sphere".
Cylinder requires divisionsU; sphere requires divisionsU and divisionsV. """
function gen_mesh(outfile, geom, divisionsU=0, divisionsV=0)
    if geom == "cube"
        mesh = cube_mesh()
    elseif geom == "cylinder"
        mesh = cylinder_mesh(divisionsU)
    elseif geom == "sphere"
        mesh = sphere_mesh(divisionsU, divisionsV)
    end
    write_obj(outfile, mesh)
end


""" est_normals(outfile, infile)
Estimate normals of the mesh stored in infile, saving the result in outfile."""
function est_normals(outfile, infile)
    input_mesh = read_obj(infile)
    mesh = estimate_normals(input_mesh)
    write_obj(outfile, mesh)
end


""" cube_mesh()
Return a new OBJMesh representing a 2x2x2 cube centered at the origin and
axis-aligned. """
function cube_mesh()
    positions = []
    uvs = []
    normals = []
    triangles = []
    # key to comments:
    # L/R = x = right/left
    # B/T = y = top/bottom
    # C/F = z = close/far
    push!(positions, Vec3( 1, -1, -1)) # 1 RBC
    push!(positions, Vec3( 1, -1,  1)) # 2 RBF
    push!(positions, Vec3(-1, -1,  1)) # 3 LBF
    push!(positions, Vec3(-1, -1, -1)) # 4 LBC
    push!(positions, Vec3( 1,  1, -1)) # 5 RTC
    push!(positions, Vec3( 1,  1,  1)) # 6 RTF
    push!(positions, Vec3(-1,  1,  1)) # 7 LTF
    push!(positions, Vec3(-1,  1, -1)) # 8 LTC

    # texture coordinates:
    push!(uvs, Vec2(1, 1)) # TR
    push!(uvs, Vec2(0, 1)) # TL
    push!(uvs, Vec2(0, 0)) # BL
    push!(uvs, Vec2(1, 0)) # BR

    # normals:
    push!(normals, Vec3( 1, 0, 0)) # R
    push!(normals, Vec3(-1, 0, 0)) # L
    push!(normals, Vec3( 0, 1, 0)) # U
    push!(normals, Vec3( 0,-1, 0)) # D
    push!(normals, Vec3( 0, 0, 1)) # C
    push!(normals, Vec3( 0, 0,-1)) # F

    # 8 faces, 2 triangles each
    push!(triangles, OBJTriangle([1,2,3], [1,2,3], [4,4,4])) # bottom face 1
    push!(triangles, OBJTriangle([1,3,4], [1,3,4], [4,4,4])) # bottom face 2
    push!(triangles, OBJTriangle([1,5,6], [4,1,2], [1,1,1])) # right face 1
    push!(triangles, OBJTriangle([1,6,2], [4,2,3], [1,1,1])) # right face 2
    push!(triangles, OBJTriangle([2,6,7], [4,1,2], [5,5,5])) # far face 1
    push!(triangles, OBJTriangle([2,7,3], [4,2,3], [5,5,5])) # far face 2
    push!(triangles, OBJTriangle([3,7,8], [2,3,4], [2,2,2])) # left face 1
    push!(triangles, OBJTriangle([3,8,4], [2,4,1], [2,2,2])) # left face 2
    push!(triangles, OBJTriangle([4,8,5], [2,3,4], [6,6,6])) # far face 1
    push!(triangles, OBJTriangle([4,5,1], [2,4,1], [6,6,6])) # far face 2
    push!(triangles, OBJTriangle([5,8,7], [1,2,3], [3,3,3])) # top face 1
    push!(triangles, OBJTriangle([5,7,6], [1,3,4], [3,3,3])) # top face 2

    # julia automatically returns the last value in the function:
    OBJMesh(positions, uvs, normals, triangles)

end


""" cylinder_mesh(n)
Return a new OBJMesh object approximation of a cylinder with radius 1 and
height 2, centered at the origin. The logitudinal axis is aligned with y, and
it is tesselated with n divisions arranged radially around the outer surface.
The ends of the cylinder are disc-shaped caps parallel to the xz plane. 
"""
function cylinder_mesh(divisionsU)
    positions = []
    normals = []
    uvs = []
    triangles = []

    theta = (2*pi)/divisionsU
    push!(positions, Vec3(0,1,0)) #Top center
    push!(positions, Vec3(0,-1,0)) #Bottom center
    push!(positions, Vec3(1,1,0))
    push!(positions, Vec3(1,-1,0))
    
    push!(normals, Vec3(0,1,0))
    push!(normals, Vec3(0,-1,0))
    x = 1 
    z = 0
    for i in 1:2:2*divisionsU
        rot_x = x*cos(theta) - z*sin(theta)
        rot_z = z*cos(theta) + x*sin(theta)
        normal_index = 3+div(i,2)
        push!(normals, normalize(positions[i+2] - positions[1])) 
        
        if i != 2*divisionsU-1 
            push!(positions, Vec3(rot_x, 1, rot_z))
            push!(positions, Vec3(rot_x, -1, rot_z))

            push!(triangles, OBJTriangle([1, i+4, i+2], [], [1,1,1])) #Top cap tris
            #Shell tris
            push!(triangles, OBJTriangle([i+2, i+4, i+5], [], [normal_index, (normal_index+1), (normal_index+1)])) 
            push!(triangles, OBJTriangle([i+5, i+3, i+2], [], [(normal_index+1), normal_index, normal_index]))
            push!(triangles, OBJTriangle([2, i+3, i+5], [], [2,2,2])) #Bottom cap tris
        else
            #Last "wedge" tris
            push!(triangles, OBJTriangle([1, 3, i+2], [], [1,1,1])) #Top
            #Shell tris
            push!(triangles, OBJTriangle([i+2, 3, 4], [], [normal_index, 3, 3])) 
            push!(triangles, OBJTriangle([i+2, 4, i+3], [], [normal_index, 3, normal_index]))
            push!(triangles, OBJTriangle([2, i+3, 4], [], [2,2,2])) #Bottom
        end
        x = rot_x
        z = rot_z
    end
    return OBJMesh(positions, uvs, normals, triangles)
end


""" sphere_mesh(n, m)
Create a Latitude-Longitude-tesselated approximation of a sphere with radius 1
centered at the origin. There are n divisions around the equator and m
divisions from pole to pole along each line of longitude. The North pole is at
(0,1,0), the South pole at (0,-1,0), and points on the Greenwich meridian are
in the x = 0 plane with z > 0. The u texture coordinate depends on longitude,
with u=0 at 180 degrees West and u=1 at 180 degrees East. The v texture
coordinate varies with latitude with v=0 at the South pole and v=1 at the North
pole. Normals should be normal to the ideal sphere surface. See the assignment
for a diagram and further details. """
function sphere_mesh(n, m)
    # TODO - feel free to drop in your A1 solution code here
end

"""
    estimate_normals(mesh::OBJMesh)
Estimates normals for the given mesh. Overwrites any existing normals and returns a new OBJMesh object.
"""
function estimate_normals(mesh::OBJMesh)
    mesh.normals = []
    for v in 1:length(mesh.positions) #Zero out each mesh normal
        push!(mesh.normals, Vec3(0,0,0))
    end
    for t in mesh.triangles #For each triangle t
        triangle_normal = get_triangle_normal(mesh.positions[t.positions[1]], mesh.positions[t.positions[2]], mesh.positions[t.positions[3]])
        for tv in t.positions #For each vertex tv around triangle t
            mesh.normals[tv] += triangle_normal
            push!(t.normals, tv) 
        end
    end
    for v in 1:length(mesh.positions) 
        mesh.normals[v] = normalize(mesh.normals[v])
    end
    return mesh
end

end # module MeshGen


