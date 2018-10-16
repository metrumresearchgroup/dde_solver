- [Introduction](#org4561907)
- [Build the static library](#org8b329e9)
- [Test](#org10f6b62)
- [C binding](#orgbc76bb0)
- [C++ binding](#org89cd3b6)
  - [Example](#org0b272ec)



<a id="org4561907"></a>

# Introduction

This repo provides a C/C++ binding for the DDE solver originally written by L. Shampine and S. Thompson: <http://www.radford.edu/~thompson/ffddes/> and its update <https://github.com/WarrenWeckesser/dde_solver>


<a id="org8b329e9"></a>

# Build the static library

To make `lib/libdde_solver.a`

```bash
make all
```


<a id="org10f6b62"></a>

# Test

In `tests`

```bash
make all
```

C binding tests: `sol1_test`, `sol2_test`, `sol3_test`, `sol4_test`. C++ binding tests: `dde_solver_test`.


<a id="orgbc76bb0"></a>

# C binding

`include/dde_solver.h` contains functions declarations of C bindings for Fortran implementations `DKL_1`, `DKL_2`, `DKL_3`, and `DKL_4`.


<a id="org89cd3b6"></a>

# C++ binding

`include/dde_solver.hpp` contains a class template `DdeIntegrator` as the of C++ binding interface.


<a id="org0b272ec"></a>

## Example

```c++
#include <dde_solver.hpp>
#include <dde_options.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <math.h>

using dde_solver_cc::DdeUserOption;
using dde_solver_cc::DdeSol;
using dde_solver_cc::DdeIntegrator;

static constexpr double R=3.5, M=19.0;

/*
 * Functor for RHS of DDEs
 */
struct FnDdes {
  std::vector<double> operator()(const double& t,
                                 const std::vector<double>& y,
                                 const std::vector<double>& z) const {
    return { R * y[0] * (1.0 - z[0] / M) };
  }
};

/*
 * Functor for history
 */
struct FnHistory {
  std::vector<double> operator()(const double& t) const {
    if (t == 0.0) {
      return {19.001};
    } else {
      return {19.0};
    }
  }
};

/*
 * Functor for events
 */
struct FnEvent
{
  std::vector<double> operator()(const double& t,
                                 const std::vector<double>& y,
                                 const std::vector<double>& dydt,
                                 const std::vector<double>& z) const {
    return { dydt[0], dydt[0] };
  }
};

FnDdes ddes;
FnHistory his;
FnEvent ef;

std::vector<double> delay{ 0.74 }; // a single delay

/*
 * user options. @c nullptr_t corresponding to optional
 * Fortran arguments that are not present. Here they are 
 * @c CHANGE_FCN, @c OUT_FCN, and @c user_trim_get.
 */
using Dde = DdeUserOption<FnDdes, double, FnHistory,
                          FnEvent, nullptr_t, nullptr_t, nullptr_t>;
Dde o1(ddes, delay, his, ef, nullptr, nullptr, nullptr);

// user options
o1.nvar.resize(3);
o1.nvar[0] = 1;               // n: size of DDE system
o1.nvar[1] = 1;               // nlags: number of delays
o1.nvar[2] = 2;               // nef: number of events

o1.tspan[0] = 0.0;            // time begin
o1.tspan[1] = 40.0;           // time end

o1.jumps.push_back(0.0);
o1.direction.push_back(1);
o1.direction.push_back(-1);

o1.rerr[0] = 1.E-4;           // relative error
o1.aerr[0] = 1.E-7;           // absolute error

o1.opts.interpolation = true;   // allow interpolation

// solution is in @c DdeSol.
DdeSol sol = DdeIntegrator<Dde>(o1)();
```

Where the solution structure `DdeSol` is of form

```c++
struct DdeSol {
  const int npts;                      // nb. of points
  const int flag;                      // nb. of delays
  const int ne;                        // nb. of events
  std::vector<double> t;               // time at points
  std::vector<double> y;               // solution at points
  std::vector<double> te;              // time at events
  std::vector<double> ye;              // solution at events
  std::vector<double> yoft;
  std::vector<int> stats;
  std::vector<int> ie;
  std::vector<int> ipoint;
  bool shift;
  double tshift;
}
```

For the usage of each input/output variable, see <https://www.radford.edu/~thompson/ffddes/>.
