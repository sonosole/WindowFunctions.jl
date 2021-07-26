export tukey
"""
    tukey(N::Int, r::AbstractArray) -> AbstractArray

Tukey window is a rectangular window with the first and
last r/2 percent of the samples equal to parts of a cosine.
"""
function tukey(N::Int, r::AbstractFloat=0.5)
    n = (0:N-1)/(N-1)
    w = zeros(N)
    for i = 1:N
        x = n[i]
        if 0 <= x <= r/2
            w[i] = 0.5 + 0.5*cos(2*pi/r*(x-r/2))
        elseif r/2 <= x <= 1-r/2
            w[i] = 1.0
        else
            w[i] = 0.5 + 0.5*cos(2*pi/r*(x-1+r/2))
        end
    end
    return w
end


export gauss
"""
    gauss(N::Int, a::AbstractFloat=0.5)

a is proportional to the standard deviation
"""
function gauss(N::Int, a::AbstractFloat=0.4)
    n = 0:N-1
    m = (N-1)/2
    return @. exp(-((n - m)/m/a)^2/2)
end
