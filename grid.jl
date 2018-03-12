include("wires.jl")

struct Grid{T} <:Any where {T<:AbstractFloat}
    girder_length::Array{T, 1}
    sample_points::Array{T, 2} # position (x_m,y_m,z_m)
    points::Int64 # total sample points M
    xyzs::Array{Array{T, 1}, 1} # grid x,y,z
    function Grid(g::T, srl::Array{T, 1}, samples::Array{I, 1}, center::Array{T, 1}) where {T, I<:Integer}
        # srl: sample region length
        # samples: including both end
        bounds = zeros(3, 2)

        bounds[:, 1] = center-srl/2
        bounds[:, 2] = center+srl/2
        # sample_space = zeros(T, 3)
        xyzs = [] 
        
        for i in 1:3
            append!(xyzs, [collect(linspace(bounds[i, 1], bounds[i, 2], samples[i]))])
        end

        sample_points = [[x, y, z] for x in xyzs[1] for y in xyzs[2] for z in xyzs[3]]
        @show size(sample_points, 1)
        new{T}([g,g,g], hcat(sample_points...), size(sample_points, 1), xyzs)
    end
end

G = Grid(1., [0., .4, .4], [1, 21, 21], 0.*ones(3))

# center: x,y = 0.5, 0.0<z<0.5/0.5<z'<1
# normal = [0, 0, 1]
# radius: r<0.5

hf_wires = 3
rzs =[0.34, 0.41, 0.5]
# rzs =[0.48, 0.49, 0.50]

# hf_wires = 5
# rzs =[0.4, 0.42, 0.44, 0.46, 0.48]

# B_j_m_2n
sample_B_wire = zeros(3, G.points, hf_wires*2)
# B^0_3_m
B_object = repeat([0., 0., 1.], outer=G.points)

#　cicular wire, radiusr = r, center = (0.5,0.5,z), parallel to Oxy
function CW_fromrz(r, z)
    return CircularWire([0., 0., z], [0., 0., 1.], r)
end

function init(G::Grid, hf_wires::Integer, rzs::Array, 
            sample_B_wire::Array{Float64, 3})
    for i in 1:hf_wires
        cw = CW_fromrz(rzs[i], G.girder_length[3]/2)
        for j in 1:G.points
            sample_B_wire[:, j, i] .= calcB(cw, G.sample_points[:, j])
        end
    end
    for i in 1:hf_wires
        cw = CW_fromrz(rzs[i], -G.girder_length[3]/2)
        for j in 1:G.points
            sample_B_wire[:, j, i+hf_wires] .= calcB(cw, G.sample_points[:, j])
        end
    end
end



function evalB(G, B_object, sample_B_wire)
    Bji = zeros(3*G.points, size(sample_B_wire, 3))
    B_object = collect(Iterators.flatten(B_object))
    for j in 1:G.points
        for i in 1:size(sample_B_wire, 3)
            Bji[(j-1)*3+1:(j-1)*3+3, i] = sample_B_wire[:, j, i]
        end
    end

    I = (Bji.'*Bji)\(Bji.'*B_object)    
    # I = 1/(4*π/(0.5*(5/4)^(3/2)))*ones(2*hf_wires) # helmholtz coil
    B_closest = Bji*I
    B_diff =  reshape((B_object - B_closest), 3, G.points)
    B_diff_mag = [norm(B_diff[:, i]) for i in 1:G.points] # average |\Delta B|
    return norm(B_diff_mag)/sqrt(size(B_diff_mag, 1)), B_diff, I
    # rms, ..., ...
end

init(G, hf_wires, rzs, sample_B_wire)#, CW_list, gCW_list)

diffBnorm, B_diff, I = evalB(G, B_object, sample_B_wire)
@show 0, diffBnorm, maximum(B_diff), I
# include("wires.jl")

function update(rzs, sample_B_wire, G, hf_wires, r_step=0.01)
    i, j, k, l = rand(4)
    i, j, k, l = floor.(Int, [i*size(rzs, 1)+1, j*3+1, k*15+1, l*2])

    Δr = r_step*(k-8)

    r = rzs[i]+Δr
    r_max = G.girder_length[1]/2
    if (r>r_max) | (r<0) | (findfirst(x->r≈x, rzs)!=0) # better use approx compare
        return update(rzs, sample_B_wire, G, hf_wires, r_step)
    end
    rzs[i] = r
    cw = CW_fromrz(rzs[i], G.girder_length[3]/2)
    for j in 1:G.points
        sample_B_wire[:, j, i] .= calcB(cw, G.sample_points[:, j])
    end

    cw = CW_fromrz(rzs[i], -G.girder_length[3]/2)
    for j in 1:G.points
        sample_B_wire[:, j, i+hf_wires] .= calcB(cw, G.sample_points[:, j])
    end
    

    sorted_arg = collect(1:length(rzs))
    sort!(sorted_arg, by=x->rzs[x])
    new_rzs = similar(rzs)
    new_sample_B_wire = similar(sample_B_wire)
    for i in 1:length(rzs)
        new_rzs[i] = rzs[sorted_arg[i]]
        new_sample_B_wire[:, :, sorted_arg[i]] = sample_B_wire[:, :, sorted_arg[i]]
        new_sample_B_wire[:, :, sorted_arg[i]+hf_wires] = sample_B_wire[:, :, sorted_arg[i]+hf_wires]
    end
    return new_rzs, new_sample_B_wire
end

n = 0
last_update = n
T = diffBnorm * 2
object = sum(I.^2)*1e-5
min_state = [object,diffBnorm, rzs, sample_B_wire, I]
# diffBnorm, rzs, sample_B_wire, I = min_state
Tscale = 0.1
Trate = 0.96
println(n, rzs)
# while n<4000
#     new_rzs, new_sample_B_wire = update(rzs, sample_B_wire, G, hf_wires)
#     valid, (new_diffBnorm, new_B_diff, new_I) = try
#         true, evalB(G, B_object, sample_B_wire)
#     catch Base.LinAlg.SingularException(6)
#         @show new_rzs
#         @show Base.LinAlg.SingularException(6)
#         false, (0,0,0)
#     end
#     if !valid
#         continue
#     end
#     new_object = new_diffBnorm + sum(abs.(new_I.^2))*1e-4
#     Δdiff = new_object - object
    
#     if (Δdiff < 0) | (rand() < exp((-Δdiff)/(Tscale*T)))
#         last_update = n
#         rzs, sample_B_wire = new_rzs, new_sample_B_wire
#         diffBnorm = new_diffBnorm
#         object = new_object
#         I = new_I
#         T = T*0.8 + (new_object) * 0.2
    
#         @show n,new_object, diffBnorm, maximum(new_B_diff), sum(abs.(new_I.^2))*1e-4
#         @show new_rzs
#         if Δdiff > 0 
#             @show Δdiff/object, exp((-Δdiff)/(Tscale*T))
#         end
#     end
    
#     if object < min_state[1]
#         min_state = [object, diffBnorm, rzs, sample_B_wire, I]
#     end
#     n += 1
#     if n%50 == 0
#         @show n
#     end
# end

using JLD
save("min_state_3_z0.jld", "min_state", min_state, "G", G, "I", I)

