export tukey
"""
    tukey(N::Int, r::AbstractFloat=0.5; dtype::DataType=Float32) -> AbstractArray

Tukey window is a rectangular window with the first and
last r/2 percent of the samples equal to parts of a cosine.
"""
function tukey(N::Int, r::AbstractFloat=0.5; dtype::DataType=Float32)
    @assert r ≤ 1.0
    R = r/2
    n = (0:N-1)/(N-1)
    w = zeros(dtype, N)
    for i = 1:N
        x = n[i]
        if 0 ≤ x ≤ R
            w[i] = 0.5 + 0.5*cos(2*pi/r*(x-R))
        elseif R ≤ x ≤ 1-R
            w[i] = 1.0
        else
            w[i] = 0.5 + 0.5*cos(2*pi/r*(x-1+R))
        end
    end
    return w
end


export gauss
"""
    gauss(N::Int, a::AbstractFloat=0.5; dtype::DataType=Float32)

a is proportional to the standard deviation
"""
function gauss(N::Int, a::AbstractFloat=0.4; dtype::DataType=Float32)
    m = (N-1)/2
    w = zeros(dtype, N)
    for i in 0:N-1
        w[i+1] = exp(-((i - m)/m/a)^2/2)
    end
    return w
end


export Γ, I₀, Iᵦ, kaiser
"""
    Γ(n::Integer)
Γ(n) = (n-1)!
"""
function Γ(n::Integer)
    x = BigInt(1)
    for i = 2:n-1
        x *= i
    end
    return x
end


"""
    I₀(x::Real, N::Integer) -> s::Real
Zero order modified Bessel function of the first kind.
    I₀(x) = ∑ᵢ x²ⁱ / Γ(i+1)², i = 0:1:Inf
"""
function I₀(x::Real, N::Integer)
    s  = 0
    x /= 2
    for n in 0:N
        s += x^(2*n) / Γ(n+1)^2
    end
    return s
end


"""
    Iᵦ(x::Real, β::Integer, N::Integer) -> s::Real
Iᵦ is the β order modified Bessel function of the first kind.
    Iᵦ(x) = ∑ᵢ x²ⁱ⁺ᵝ / (Γ(i+1) * Γ(i+β+1)), i = 0:1:Inf
"""
function Iᵦ(x::Real, β::Integer, N::Integer)
    s  = 0
    x /= 2
    for n in 0:N
        s += x^(2*n + β) / (Γ(n+1) * Γ(n+β+1))
    end
    return s
end


"""
    kaiser(N::Integer, α::Real, M::Integer; dtype::DataType=Float32) -> w::Array
The Kaiser window is an approximation to the prolate spheroidal window, for\n
which the ratio of the mainlobe energy to the sidelobe energy is maximized.
α is any real number
"""
function kaiser(N::Integer, α::Real, M::Integer; dtype::DataType=Float32)
    πα = π * α
    w  = zeros(dtype, N)
    D  = I₀(πα, M)
    for n = 0:N-1
        w[n+1] = I₀(πα * sqrt(1 - (2n/(N-1) - 1)^2), M) / D
    end
    return w
end

