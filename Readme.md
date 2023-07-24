This repository contains implementations of the semi interleaved ladder to measure its performances and some other tools :

* **intel_measurement_stuff.c** : tools to measure cycle count and the number of instructions
* **measure-ladder-mont-semi.c** : source code of the semi interleaved ladder and the Montgomery ladder using GnuMP low level functions,
* **measure-ladder-par.c** : source code to get the performances of the semi interleaved ladder using 2 cores,
* **measure-ladder-par3.c** : source code to get the performances of the semi interleaved ladder using 3 cores.
* **measure.sh** : bash script to run before any test. It allows to count the instructions of a running program and it disables the Turbo boost feature so as to obtain stable measures.

**Test Platform** : All the codes have been tested on a Intel(R) Core(TM) i9-11900KF @ 3.50GHz, gcc 11.3.0, Ubuntu 22.04

**Prerequisites**

To run the tests you must have :
* the GnuMP library
* the OpenMP package (for 2 cores and 3 cores version of the semi interleaved ladder)
* the msr-tools

**How to run ?**

First configure the msr-tools and disable the turbo-boost feature. In a shell, run:
```console
sudo bash measure.sh
```

Next as a normal user, execute :

```console
gcc -03 measure-ladder-mont-semi.c -O measure-ladder-mont-semi -lgmp
./measure-ladder-mont-semi size_in_bits
```
to get a comparison of the performances of the semi interleaved ladder and the classical Montgomery ladder. The value *size_in_bits* is the size of the values used to compute a modular exponentiation using the semi interleaved ladder and the Montgomery ladder .

To use the parallel version, execute (2 core version) :

```console
gcc -03 measure-ladder-par.c -O measure-ladder-par -fopenmp -lgmp
./measure-ladder-par size_in_bits
```

or (3 core version) :

```console
gcc -03 measure-ladder-par3.c -O measure-ladder-par -fopenmp -lgmp
./measure-ladder-par3 size_in_bits
```


