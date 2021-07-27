"""
In both digital filter design and spectral estimation, the choice of a windowing
function can play an important role in determining the quality of overall results.
The main role of the window is to damp out the effects of the Gibbs phenomenon that
results from truncation of an infinite series.

# Available Window Functions

    | Function       | Window Type                                    |
    | -------------- | ---------------------------------------------- |
    | barthann       | Modified Bartlett-Hann                         |
    | bartlett       | Bartlett                                       |
    | blackman       | Blackman                                       |
    | blackmanharris | Minimum four-term Blackman-Harris              |
    | bohman         | Bohman                                         |
    | cheb           | Chebyshev                                      |
    | flattop        | Flat top weighted                              |
    | hamming        | Hamming                                        |
    | hanning        | Hann (Hanning)                                 |
    | nuttall        | Nuttall-defined minimum 4-term Blackman-Harris |
    | parzen         | Parzen (de la Vallée Poussin)                  |
    | rectangular    | Rectangular                                    |
    | triangular     | Triangular                                     |
    | tukey          | Tukey (tapered cosine)                         |
    | taylor         | Taylor                                         |
    | gauss          | Gaussian                                       |
    | kaiser         | Kaiser                                         |
    | enbw           | Equivalent noise bandwidth                     |

"""
module WindowFunctions

include("./1-unary.jl")
include("./2-binary.jl")

end # WindowFunctions
