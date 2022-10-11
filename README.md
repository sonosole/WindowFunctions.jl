
# WindowFunctions

## Motivation

In both digital filter design and spectral estimation, the choice of a windowing
function can play an important role in determining the quality of overall results.
The main role of the window is to damp out the effects of the Gibbs phenomenon that
resulted from truncation of an infinite series.

## Methods

| Function       | Window Type                                    |
| -------------- | ---------------------------------------------- |
| barthann       | Modified Bartlett-Hann                         |
| bartlett       | Bartlett                                       |
| blackman       | Blackman                                       |
| blackmanharris | Minimum four-term Blackman-Harris              |
| bohman         | Bohman                                         |
| ~~cheb~~           | Chebyshev                                      |
| flattop        | Flat top weighted                              |
| hamming        | Hamming                                        |
| hanning        | Hann (Hanning)                                 |
| nuttall        | Nuttall-defined minimum 4-term Blackman-Harris |
| parzen         | Parzen (de la Vallée Poussin)                  |
| rectangular    | Rectangular                                    |
| triangular     | Triangular                                     |
| tukey          | Tukey (tapered cosine)                         |
| ~~taylor~~         | Taylor                                         |
| gauss          | Gaussian                                       |
| kaiser         | Kaiser                                         |
| ~~enbw~~           | Equivalent noise bandwidth                     |


## Installation
To install the released stable version, run
```julia
using Pkg
Pkg.add("WindowFunctions")
```

## Example

```julia
using Plots
using WindowFunctions

kaiser5(N; dtype::DataType=Float32) = kaiser(N, 5.0, 10, dtype=dtype)

winlist = [barthann,
           bartlett,
           blackman,
           blackmanharris,
           bohman,
           flattop,
           hamming,
           hanning,
           nuttall,
           parzen,
           rectangular,
           triangular,
           gauss,
           tukey,
           kaiser5]

for win in winlist
    plot(win(256, dtype=Float32),
         label=string(win),
         fillrange = 0,
         fillalpha = 0.3,
         show=true)
    ylims!(-0.13,1.18)
    xlims!(1-20,256+20)
    sleep(1)
end
```
