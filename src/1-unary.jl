export barthann
"""
    barthann(N::Int; dtype::DataType=Float32) -> AbstractArray
    Sidelobe peak attenuation = 35.9dB
barthann window has a mainlobe at the origin and asymptotically decaying
sidelobes on both sides. It is a linear combination of weighted Bartlett
and Hanning windows with near sidelobes lower than both Bartlett and Hanning
and with far sidelobes lower than both Bartlett and Hamming windows. The
mainlobe width of the modified Bartlett-Hann window is not increased relative
to either Bartlett or Hann window mainlobes.
"""
function barthann(N::Int; dtype::DataType=Float32)
    w = Array{dtype}(undef, N)
    for n = 1:N
        x = (n-1)/(N-1)
        w[n] = 0.62 - 0.48*abs(x-0.5) + 0.38*cos(2π*x - π)
    end
    return w
end


export bartlett
"""
    bartlett(N::Int; dtype::DataType=Float32) -> AbstractArray
    Sidelobe peak attenuation = 26.5dB
Bartlett window is very similar to a triangular window as returned by the triang
function. However, the Bartlett window always has zeros at the first and last samples,
while the triangular window is nonzero at those points.
"""
function bartlett(L::Int; dtype::DataType=Float32)
    L==1 && return dtype[1.0]
    N = L - 1
    w = zeros(dtype, L)
    for i = 1:L
        n = i - 1
        w[i] = 0<=n<=N/2 ? 2n/N : 2-2n/N
    end
    return w
end


export blackman
"""
    blackman(N::Int; dtype::DataType=Float32) -> AbstractArray
    Sidelobe peak attenuation = 58.1dB
It is recommended to use blackman window to detect two
signals with similar frequency but different amplitude.
"""
function blackman(N::Int; dtype::DataType=Float32)
    w = Array{dtype}(undef, N)
    k = 2*π / (N-1)
    for n = 1:N
        x = k * (n-1)
        w[n] = 0.42 - 0.5 * cos(x) + 0.08 * cos(2x)
    end
	return w
end


export blackmanharris
"""
    blackmanharris(N::Int; dtype::DataType=Float32) -> AbstractArray
    Sidelobe peak attenuation = 92dB
Minimum four-term Blackman-Harris window.
"""
function blackmanharris(N::Int; dtype::DataType=Float32)
    a₀ = 0.35875
    a₁ = 0.48829
    a₂ = 0.14128
    a₃ = 0.01168
    w = Array{dtype}(undef, N)
    k = 2*π / (N-1)
    for n = 1:N
        x = k * (n-1)
        w[n] = a₀ - a₁*cos(x) + a₂*cos(2x) - a₃*cos(3x)
    end
    return w
end


export bohman
"""
    bohman(N::Int; dtype::DataType=Float32) -> AbstractArray
    Sidelobe peak attenuation = 46dB
Bohman window is the convolution of two half-duration cosine lobes.
In the time domain, it is the product of a triangular window and a
single cycle of a cosine with a term added to set the first derivative
to zero at the boundary. Bohman windows fall off as 1/ω⁴
"""
function bohman(N::Int; dtype::DataType=Float32)
    w = Array{dtype}(undef, N)
    for i = 1:N
        n = (i-1)/(N-1)
        n = abs(2n - 1)
        w[i] = (1-n)*cos(π*n) + 1/π*sin(π*n)
    end
    return w
end


export flattop
"""
    flattop(N::Int; dtype::DataType=Float32) -> AbstractArray
    Sidelobe peak attenuation = 93.6dB
flat-top window is summations of cosines with very small passband fluctuations.
"""
function flattop(N::Int; dtype::DataType=Float32)
    a₀ = 0.215578950
    a₁ = 0.416631580
    a₂ = 0.277263158
    a₃ = 0.083578947
    a₄ = 0.006947368
    w = Array{dtype}(undef, N)
    k = 2*π / (N-1)
    for n = 1:N
        x = k * (n - 1)
        w[n] = a₀ - a₁*cos(x) + a₂*cos(2x) - a₃*cos(3x) + a₄*cos(4x)
    end
    return w
end


export hamming
"""
    hamming(N::Int; dtype::DataType=Float32) -> AbstractArray
    Sidelobe peak attenuation = 43.2dB
Similar to hanning but has lower Sidelobe.
"""
function hamming(N::Int; dtype::DataType=Float32)
    w = Array{dtype}(undef, N)
    k = 2*π / (N-1)
    for n = 1:N
        x = n - 1
        w[n] = 0.54 - 0.46 * cos(k * x)
    end
	return w
end


export hanning
"""
    hanning(N::Int; dtype::DataType=Float32) -> AbstractArray
    Sidelobe peak attenuation = 31.5dB

If the signal is random or unknown, or has multiple frequency components, and
the test focuses on frequency points rather than energy, it is recommended to
select hanning window.
"""
function hanning(N::Int; dtype::DataType=Float32)
    w = Array{dtype}(undef, N)
    k = 2*π / (N-1)
    for n = 1:N
        x = n - 1
        w[n] = 0.5 - 0.5 * cos(k * x)
    end
	return w
end


export nuttall
"""
    nuttall(N::Int; dtype::DataType=Float32) -> AbstractArray
    Sidelobe peak attenuation = 93.8dB
The window is minimum in the sense that its maximum sidelobes are minimized.
The coefficients for this window differ from the Blackman-Harris window coefficients
computed with blackmanharris and produce slightly lower sidelobes.
"""
function nuttall(N::Int; dtype::DataType=Float32)
    a₀ = 0.3635819
    a₁ = 0.4891775
    a₂ = 0.1365995
    a₃ = 0.0106411
    w = Array{dtype}(undef, N)
    k = 2*π / (N-1)
    for n = 1:N
        x = k * (n - 1)
        w[n] = a₀ - a₁*cos(x) + a₂*cos(2x) - a₃*cos(3x)
    end
    return w
end


export parzen
"""
    parzen(N::Int; dtype::DataType=Float32) -> AbstractArray
    Sidelobe peak attenuation = 53.1dB
Parzen windows are piecewise-cubic approximations of Gaussian windows.
Parzen window sidelobes fall off as 1/ω⁴
"""
function parzen(N::Int; dtype::DataType=Float32)
    w = Array{dtype}(undef, N)
    for i = 1:N
        m = i - 1
        n = abs(m - (N-1)/2)
        w[i] = n<(N-1)/4 ? 1 - 6(2n/N)^2 + 6(2n/N)^3 : 2(1-2n/N)^3
    end
    return w
end


export rectangular
"""
    rectangular(N::Int; dtype::DataType=Float32) -> AbstractArray
    Sidelobe peak attenuation = 13.3dB
"""
function rectangular(N::Int; dtype::DataType=Float32)
    return ones(dtype, N)
end


export triangular
"""
    triangular(N::Int; dtype::DataType=Float32) -> AbstractArray
    Sidelobe peak attenuation = 26.5dB
Compared with a rectangular window, the mainlobe is about twice wide, but the
sidelobes are smaller and have no negative sidelobes.
"""
function triangular(N::Int; dtype::DataType=Float32)
    N%2==0 && return triangulareven(N, dtype)
    N%2==1 && return triangularodd(N, dtype)
end

function triangularodd(N::Int, dtype::DataType=Float32)
    w = Array{dtype}(undef, N)
    for n = 1:N
        w[n] = 1 ≤ n ≤ (N+1)/2 ? 2n/(N+1) : 2-2n/(N+1)
    end
    return w
end

function triangulareven(N::Int, dtype::DataType=Float32)
    w = Array{dtype}(undef, N)
    for n = 1:N
        w[n] = 1 ≤ n ≤ N/2 ? (2n-1)/N : 2-(2n-1)/N
    end
    return w
end
