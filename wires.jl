using FastGaussQuadrature

"abstract type Wire. All concrete wire types should be its descendants."
abstract type Wire end


"""
    StraightWire(p1::Array{T}, p2::Array{T})

`p1`, `p2` are positions of two end of the straight wire in 3-D rectangular-coordinate system.
"""
struct StraightWire{T} <: Wire where T<:AbstractFloat
    p1::Array{T}
    p2::Array{T}
end

"""
    calcB(sw::StraightWire{T}, p0::Array{T, 1})

Calculate the dimentionless unit magnetic field at position `p0` of the straight wire `sw`. Coefficient μ_0*I/4π is neglated. 
"""
function calcB(sw::StraightWire{T}, p0::Array{T, 1}) where {T}
    p1, p2 = sw.p1, sw.p2
    p2_p1 = p2-p1
    tp = (vecdot(p1, p1) - vecdot(p0, p1) - vecdot(p1, p2) + vecdot(p0, p2))/(vecdot(p2_p1, p2_p1))
    pp = p1 + tp*p2_p1
    rp = p0 - pp
    
    θ1 = atan2(sign(vecdot(p1-pp, p2-p1))*(norm(p1-pp)), norm(rp))
    θ2 = atan2(sign(vecdot(p2-pp, p2-p1))*(norm(p2-pp)), norm(rp))

    n = cross(p2 - p1, p0 - p2)
    n = n/norm(n)

    
    return n*(sin(θ2)-sin(θ1))/norm(rp) # *μ_0*I/4π
end

"""
    parametricEq(sw::StraightWire)

Return the parametric equation f(t) of StraightWire `sw`, t ∈ [0, 1]
"""
function parametricEq(sw::StraightWire)
    return t->sw.p1 + t*(sw.p2-sw.p1)
end

"""
    CircularWire(center::Array{T, 1}
        normal::Array{T, 1}
        radius::T)

`center` is the center of the circular wire in 3-D rectangular-coordinate system.
`normal` is the normal vector of the circular wire in 3-D rectangular-coordinate system.
`radius` is the radius of the circular wire.
"""
struct CircularWire{T} <: Wire where T<:AbstractFloat
    center::Array{T, 1}
    normal::Array{T, 1}
    radius::T
    curve::Function
    ft::Function
end

function CircularWire(center::Array{T, 1}, normal::Array{T, 1}, radius::T) where{T<:AbstractFloat}
    z = [0., 0., 1.]
    
    x2 = cross(z, normal)
    y2 = cross(normal, x2)

    if norm(x2) ≈ 0
        x2 = [1., 0., 0.]
        y2 = [0., 1., 0.]
    else
        x2 = x2/norm(x2)
        y2 = y2/norm(y2)
    end

    curve(t)  = radius*(cos(t)*x2 + sin(t)*y2) + center

    function ft(t, p0::Array{T, 1})
        p1 = curve(t)
        dl = -radius*sin(t)*x2 + radius*cos(t)*y2
        r = p0 - p1
        return cross(dl, r)/norm(r)^3
    end
    CircularWire{T}(center, normal/norm(normal), radius, curve, ft)
end


nodes, weights = gausslegendre(300)
"""
    calcB(cw::CircularWire, p0, nodes=nodes, weights=weights)

Special method for calculate the dimentionless unit magnetic field at position `p0` of the CircularWire `cw`. 
Coefficient μ_0*I/4π is neglated. 
Use 300 points Gauss-Legendre intergral by default.
"""
function calcB(cw::CircularWire, p0, nodes=nodes, weights=weights)
    ts = nodes*π # -\pi ~ \pi
    fts = hcat([cw.ft(t, p0) for t in ts]...)
    return π*[dot(weights, fts[i, :])for i in 1:3]
end

"""
    GeneralWire(curve::Function, t1::T, t2::T)

General wire described by parametric curve function `curve`:R->R^3 with upper and lower bounds t1 and t2,
and output a 3-D vector in rectangular-coordinate system.
"""
struct GeneralWire{T} <: Wire where T<:AbstractFloat
    curve::Function
    t1::T
    t2::T
end

"""
    calcB(gw::GeneralWire{T}, vec::Array{T, 1}, n=1000)

Calculate the dimentionless unit magnetic field at position `p0` of GeneralWire `gw`. 
Coefficient μ_0*I/4π is neglated. 
Use 1000 points Gauss-Legendre of Biot-Savart intergral by default.
"""
function calcB(gw::GeneralWire{T}, vec::Array{T, 1}, n=1000) where {T}
    # nodes, weights = gausslegendre(n)
    p1 = gw.curve(gw.t1)
    step = (gw.t2-gw.t1)/n
    t = gw.t1
    B = 0
    for i in 1:n
        t = t + step
        p2 = gw.curve(t)
        dl = p2 - p1
        p1 = p2
        R = (vec - (p2+p1)/2)
        B += cross(dl, R)/norm(R)^3
    end
    return B
end

if isinteractive()

    p0 = [1., 2., 3.]
    p1 = [0., 0., 0.]
    p2 = [4., 4., 4.]



    f(t) = p1 + t*(p2-p1)
    @elapsed calcB(CircularWire([1., 0., 0.], [1., 0., 0.], 1.), [0., 0., 0.])
    @show (calcB(StraightWire(p1, p2), p0) - calcB(GeneralWire(f, 0.0, 1.0), p0))
    @show (2*π*1.^2/(1.^2+1.^2)^(3/2) - calcB(CircularWire([1., 0., 0.], [1., 0., 0.], 1.), [0., 0., 0.])[1]) # 2.28e-7
end