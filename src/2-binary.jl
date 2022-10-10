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
