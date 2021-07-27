export barthann
"""
    barthann(N::Int) -> AbstractArray
    Sidelobe peak attenuation = 35.9dB
barthann window has a mainlobe at the origin and asymptotically decaying
sidelobes on both sides. It is a linear combination of weighted Bartlett
and Hanning windows with near sidelobes lower than both Bartlett and Hanning
and with far sidelobes lower than both Bartlett and Hamming windows. The
mainlobe width of the modified Bartlett-Hann window is not increased relative
to either Bartlett or Hann window mainlobes.
"""
function barthann(N::Int)
    n = (0:N-1)/(N-1)
    return @. 0.62 - 0.48*abs(n-0.5) + 0.38*cos(2*pi*n-pi)
end


export bartlett
"""
    bartlett(N::Int) -> AbstractArray
    Sidelobe peak attenuation = 26.5dB
Bartlett window is very similar to a triangular window as returned by the triang
function. However, the Bartlett window always has zeros at the first and last samples,
while the triangular window is nonzero at those points.
"""
function bartlett(L::Int)
    L==1 && return [1.0]
    N = L - 1
    w = zeros(L)
    for i = 1:L
        n = i - 1
        w[i] = 0<=n<=N/2 ? 2n/N : 2-2n/N
    end
    return w
end


export blackman
"""
    blackman(N::Int) -> AbstractArray
    Sidelobe peak attenuation = 58.1dB
It is recommended to use blackman window to detect two
signals with similar frequency but different amplitude.
"""
function blackman(N::Int)
    n = 2 * pi .* (0:(N-1))/(N-1)
	return @. 0.42 - 0.5 * cos(n) + 0.08 * cos(2n)
end


export blackmanharris
"""
    blackmanharris(N::Int) -> AbstractArray
    Sidelobe peak attenuation = 92dB
Minimum four-term Blackman-Harris window.
"""
function blackmanharris(N::Int)
    n = 2 * pi .* (0:N-1)/(N-1)
    a₀ = 0.35875
    a₁ = 0.48829
    a₂ = 0.14128
    a₃ = 0.01168
    return @. a₀ - a₁*cos(n) + a₂*cos(2n) - a₃*cos(3n)
end


export bohman
"""
    bohman(N::Int) -> AbstractArray
    Sidelobe peak attenuation = 46dB
Bohman window is the convolution of two half-duration cosine lobes.
In the time domain, it is the product of a triangular window and a
single cycle of a cosine with a term added to set the first derivative
to zero at the boundary. Bohman windows fall off as 1/ω⁴
"""
function bohman(N::Int)
    n = (0:N-1)/(N-1)
    n = abs.((2n) .- 1)
    return @. (1-n)*cos(pi*n) + 1/pi*sin(pi*n)
end


export flattop
"""
    flattop(N::Int) -> AbstractArray
    Sidelobe peak attenuation = 93.6dB
flat-top window is summations of cosines with very small passband fluctuations.
"""
function flattop(N::Int)
    n = 2 * pi .* (0:N-1)/(N-1)
    a₀ = 0.215578950
    a₁ = 0.416631580
    a₂ = 0.277263158
    a₃ = 0.083578947
    a₄ = 0.006947368
    return @. a₀ - a₁*cos(n) + a₂*cos(2n) - a₃*cos(3n) + a₄*cos(4n)
end


export hamming
"""
    hamming(N::Int) -> AbstractArray
    Sidelobe peak attenuation = 43.2dB
Similar to hanning but has lower Sidelobe.
"""
function hamming(N::Int)
	return @. 0.54 - 0.46 * cos( 2*pi * (0:(N-1))/(N-1) )
end


export hanning
"""
    hanning(N::Int) -> AbstractArray
    Sidelobe peak attenuation = 31.5dB

If the signal is random or unknown, or has multiple frequency components, and
the test focuses on frequency points rather than energy, it is recommended to
select hanning window.
"""
function hanning(N::Int)
	return @. 0.5 - 0.5 * cos( 2*pi * (0:(N-1))/(N-1) )
end


export nuttall
"""
    nuttall(N::Int) -> AbstractArray
    Sidelobe peak attenuation = 93.8dB
The window is minimum in the sense that its maximum sidelobes are minimized.
The coefficients for this window differ from the Blackman-Harris window coefficients
computed with blackmanharris and produce slightly lower sidelobes.
"""
function nuttall(N::Int)
    n = 2 * pi .* (0:N-1)/(N-1)
    a₀ = 0.3635819
    a₁ = 0.4891775
    a₂ = 0.1365995
    a₃ = 0.0106411
    return @. a₀ - a₁*cos(n) + a₂*cos(2n) - a₃*cos(3n)
end


export parzen
"""
    parzen(N::Int) -> AbstractArray
    Sidelobe peak attenuation = 53.1dB
Parzen windows are piecewise-cubic approximations of Gaussian windows.
Parzen window sidelobes fall off as 1/ω⁴
"""
function parzen(N::Int)
    w = zeros(N)
    for i = 1:N
        m = i - 1
        n = abs(m - (N-1)/2)
        w[i] = n<(N-1)/4 ? 1 - 6(2n/N)^2 + 6(2n/N)^3 : 2(1-2n/N)^3
    end
    return w
end


export rectangular
"""
    rectangular(N::Int) -> AbstractArray
    Sidelobe peak attenuation = 13.3dB
"""
function rectangular(N::Int)
    return ones(N)
end


export triangular
"""
    triangular(N::Int) -> AbstractArray
    Sidelobe peak attenuation = 26.5dB
Compared with a rectangular window, the mainlobe is about twice wide, but the
sidelobes are smaller and have no negative sidelobes.
"""
function triangular(N::Int)
    N%2==0 && return triangulareven(N)
    N%2==1 && return triangularodd(N)
end

function triangularodd(N::Int)
    w = zeros(N)
    for n = 1:N
        w[n] = 1<=n<=(N+1)/2 ? 2n/(N+1) : 2-2n/(N+1)
    end
    return w
end

function triangulareven(N::Int)
    w = zeros(N)
    for n = 1:N
        w[n] = 1<=n<=N/2 ? (2n-1)/N : 2-(2n-1)/N
    end
    return w
end
