# Plum

Plum is a [Monte Carlo simulation](https://en.wikipedia.org/wiki/Monte_Carlo_method) package for [polymers](https://en.wikipedia.org/wiki/Polymer). Plum contains many of the **common features** in typical Monte Carlo packges:
* Canonical ensemble Monte Carlo simulations for arbitrary linear polymer-small ion systems.
* Lennard-Jones, electrostatics and bond interactions between atoms/beads.
* Ewald summation for electrostatics calculations.
* Flexible atom/bead parameter settings using bead type and bead partial charge.
* Pivot, crankshaft, replation and center of mass translation Monte Carlo moves for the polymers.

Plum also contains the following **special features**:
* Grand-canonical ensemble simulation with configurational-bias chain insertion and deletion.
* Confined simulations between two infinite surfaces that are periodic along x and y directions.
  * End-grafted polymer brush simulations.
  * Arbitrary distribution surface interaction sites, charged and/or Lennard-Jones.
  * Uniform surface Lennard-Jones potential.
* Osmotic pressure calculations for both bulk and confined systems.

In the near future, we would also like to implement the following features in Plum:
* The routines to calculate angle and dihedral angle energies and forces.
* The routines to simulate arbitrarily branched polymers.
* The routines to handle arbitrarily grafted polymers on the confining surfaces.

## Getting Started

The following instructions explain how to compile Plum on your local machine. [Example](examples) Plum simulations are also provided to help you get started with your own project using Plum.

Note that this guide is for the Linux environment only, but it should also apply to the Mac OS terminal and Xcode.

### Prerequisites

Plum uses free C++ library Eigen and requires Eigen during compilation. Before compiling Plum, [download Eigen](https://eigen.tuxfamily.org/) and decompress it to your local directory. Plum has been tested to work with Eigen version 3.2.8.

You should also make sure that `cmake` and `gcc` are installed on your machine.

### Installation

Download the entire [src](src) directory from the Plum repository, store in `your_path_to_plum/plum`. Do

```
cd your_path_to_plum/plum
mkdir bin
```

And then change line

```
INC=-I /home/nuowang/bin/eigen-3.2.8
```

in `your_path_to_plum/plum/src/Makefile` to

```
INC=-I your_path_to_eigen/eigen_x.x.x
```

"x.x.x" represents the version of Eigen that you are using. Last

```
make
```

and you will find the binary executable `plum` under `your_path_to_plum/plum/bin`

## Runing examples

You will find 4 example simulation runs in the [examples](examples) directory. Run each example by

```
your_path_to_plum/bin/plum < run.in > run.log 
```

## Contributing

Currently, we have not set up the rules for submitting pull requests. Please directly contact [the developer](https://github.com/nuowang) if you would like to contribute.

## Authors

* **Nuo Wang** - The current developer of Plum, responsible for most of the functionalities in the current version of the code.
* **Rachel Krueger** - The inital framework of the code.

## License

Plum is licensed under the MIT License, see [LICENSE.md](LICENSE.md) for details.

## Acknowledgments

* Dr. Pengfei Zhang for his guidance, suggestions and help in code validation.
* Rachel Krueger for creating the initial framework of what became Plum.
* Dow Chemical Co. for funding the research projects for which Plum was written.
