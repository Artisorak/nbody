# N Body Simulation

This code simulates particles under self-gravity. The forces are calculated via direct summation or with a tree and multipole expansion. 

The implementation in this branch also uses adaptive timestepping. I mostly worked on the non-adaptive implementation though, so it's possible that I forgot to implement some of the improvements I made for this branch. 

Most of the code is written in C++, Python is used to make plots and animations. OpenMP is used for parallelisation. 

## Install OpenMP

https://www.geeksforgeeks.org/openmp-introduction-with-installation-guide/

## Clone repository

This will create a folder called nbody that contains everything. 

```
git clone https://github.com/Artisorak/nbody.git
```

## Make build directory

This requires cmake and make. Run it in the folder nbody. 

```
mkdir build
cd build
cmake ..
```

## Compiling and running the program

Go to the nbody/build directory and run these commands: 

```
make
./nbody
```
