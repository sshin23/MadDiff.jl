# MadDiff.jl


| **Documentation** | **Build Status** | **Coverage** |
|:-----------------:|:----------------:|:----------------:|
| [![doc](https://img.shields.io/badge/docs-dev-blue.svg)](https://sshin23.github.io/MadDiff.jl/dev) | [![build](https://github.com/sshin23/MadDiff.jl/actions/workflows/test.yml/badge.svg)](https://github.com/sshin23/MadDiff.jl/actions/workflows/test.yml) | [![codecov](https://codecov.io/gh/sshin23/MadDiff.jl/branch/main/graph/badge.svg?token=U6NMMW0IT5)](https://codecov.io/gh/sshin23/MadDiff.jl) |



**This package is under development. Use at your own risk.**

MadDiff.jl is an automatic differentiation and algebraic modeling package. MadDiff.jl constructs first and second derivative functions off-line (i.e., prior to calling the optimization solver) by applying operator overloading-based automatic differentiation on _functions_. The exact _derivative functions_ can be obtained as results. A benchmark code can be found [here](https://github.com/sshin23/MadDiff.jl/tree/main/benchmark).

## Bug reports and support
Please report issues and feature requests via the [Github issue tracker](https://github.com/sshin23/MadDiff.jl/issues). 
