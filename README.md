# multirate-mpc-cbf

This repo containes the multirate-mpc-cbf package, which is a layered multi-frequency control architecrure introduced in [1].

## Prerequisite 

The [OSQP](https://github.com/oxfordcontrol/osqp) library is required.

## Installation

This library can be installed using the following commands:

```
mkdir build && cd build
cmake ..
make install
```

## Unit Test

To check if the library is installed correctly you can test the MPC
```
reset
g++ src/mpc_test.cpp -o mpc_test.o -I /usr/local/include/osqp/ -L /usr/local/lib/osqp/ -losqp -L /usr/local/lib/mpc/ -lmpc
./mpc_test.o
```
and the CBF
```
reset
g++ src/cbf_test.cpp -o cbf_test.o -I /usr/local/include/osqp/ -L /usr/local/lib/osqp/ -losqp -L /usr/local/lib/mpc/ -lmpc
./cbf_test.o
```


## Integration in ROS

Checkout our [Segway](https://github.com/DrewSingletary/segway_sim) ROS simulator to see how to used this library.

## References

This code is based on the following article:

1. U. Rosolia and A. D. Ames, "Multi-Rate Control Design Leveraging Control Barrier Functions and Model Predictive Control Policies," in IEEE Control Systems Letters, vol. 5, no. 3, pp. 1007-1012, July 2021, [PDF](https://ieeexplore.ieee.org/document/9137248)
