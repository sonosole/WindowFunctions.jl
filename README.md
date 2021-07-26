
# WindowFunctions

## Motivation

In both digital filter design and spectral estimation, the choice of a windowing
function can play an important role in determining the quality of overall results.
The main role of the window is to damp out the effects of the Gibbs phenomenon that
results from truncation of an infinite series.

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
| ~~kaiser~~         | Kaiser                                         |
| ~~enbw~~           | Equivalent noise bandwidth                     |


## Installation
To install the released stable version, enter the REPL mode
```julia
] add WindowFunctions
```
or
```julia
using Pkg
Pkg.add("WindowFunctions")
```

## Example

```julia
using Plots
using WindowFunctions

plot(gauss(256,0.20),label="gaussian 0.20")
plot!(gauss(256,0.35),label="gaussian 0.35")
plot!(flattop(256),label="flattop")
```
