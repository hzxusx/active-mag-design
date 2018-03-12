# struct A{T}
#     x::T
#     y::T
#     f::Function
# end

# function A(x)
#     f(t) = x+t 
#     A(x, x, f)
# end

# @show A(1)

# using FastGaussQuadrature
# struct CircularWire{T} <:Any where T<:AbstractFloat
#     center::Array{T, 1}
#     normal::Array{T, 1}
#     radius::T
#     curve::Function
#     ft::Function
# end

# function CircularWire(center::Array{T, 1}, normal::Array{T, 1}, radius::T) where{T<:AbstractFloat}
#     z = [0., 0., 1.]
    
#     x2 = cross(z, normal)
#     y2 = cross(normal, x2)

#     if norm(x2) ≈ 0
#         x2 = [1., 0., 0.]
#         y2 = [0., 1., 0.]
#     else
#         x2 = x2/norm(x2)
#         y2 = y2/norm(y2)
#     end

#     curve(t)  = radius*(cos(t)*x2 + sin(t)*y2) + center

#     function ft(t, p0::Array{T, 1})
#         p1 = curve(t)
#         dl = -radius*sin(t)*x2 + radius*cos(t)*y2
#         r = p0 - p1
#         return cross(dl, r)/norm(r)^3
#     end
#     CircularWire{T}(center, normal/norm(normal), radius, curve, ft)
# end


# nodes, weights = gausslegendre(300)
# function calcB(cw::CircularWire, p0, nodes=nodes, weights=weights)
#     ts = nodes*π
#     fts = hcat([cw.ft(t, p0) for t in ts]...)
#     @show size(fts)
#     return π*[dot(weights, fts[i, :])for i in 1:3]
# end

# @show calcB(CircularWire([1., 0., 0.], [1., 0., 0.], 1.), [0., 0., 0.])
# @show 2*π*1.^2/(1.^2+1.^2)^(3/2)
# @show (2*π*1.^2/(1.^2+1.^2)^(3/2) - calcB(CircularWire([1., 0., 0.], [1., 0., 0.], 1.), [0., 0., 0.])[1]) # 2.28e-7


try 
    x=1
catch Base.DivideError
    println(x)
end

print(x)