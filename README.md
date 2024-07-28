# DCRP

## About DCRP

DCRP means detailed cable routing problem in distributed photovoltaic systems. This is a project for solving the detailed cable routing problem in distributed photovoltaic systems using a variable-depth large neighborhood search algorithm.

## Requirements 

* CMake 3.14
* Modern C++ compiler like Visual Studio 2019 or gcc 11.4.
* (optional, only supported on Linux) Gurobi 10

## Compilation 

DCRP uses CMake for building. 
To build the project, run the following command in the root directory of the project:
```
mkdir build
cd build
cmake ..
```
Then compile the code using 
```
cmake --build . --config Release
```
This would compile the DCRP library and executable file `DCRPVDLNS.exe` (on Windows) or `DCRPVDLNS` (on other platforms) at the folder `./bin/`.


## Usage 

To run the algorithm, use the following command:
```
./bin/DCRPVDLNS <input_file> <output_file> [-s <statistic file>] [-d <max depth>]
```
where 
- `<input_file>` is the input file containing the problem instance.
- `<output_file>` is the output file where the solution will be written.
- `-s <statistic file>` is an optional argument for writing the statistics of the algorithm to a file. It is defaultly disabled.
- `-d <max depth>` is an optional argument for setting the maximum search depth. Its default value is 3.

## File Formats

### Input File

The input file should be in the following format:
```
p // number of modules
q // number of rectangular barriers
m // partition size
n_1 n_2... n_m // group sizes
x_1 y_1 // position of module 1
...
x_p y_p // position of module p
xmin_1 xmax_1 ymin_1 ymax_1 // barrier 1 = [xmin_1,xmax_1]X[ymin_1,ymax_1]
...
xmin_q xmax_q ymin_q ymax_q // barrier q = [xmin_q,xmax_q]X[ymin_q,ymax_q]
```


### Output File 

The output file would be in the following format:
```
m // partition 
n_1 n_2... n_m // group sizes
x_{11} y_{11};x_{12} y_{12}...x_{1n_1} y_{1n_1} // modules in the 1st group
...
x_{m1} y_{m1};x_{m2} y_{m2}...x_{mn_m} y_{mn_m} // modules in the mth group
```

### Statistic File 

The statistic file is a csv file containing the following columns:
- `data` path of the input file 
- `depth` maximum search depth
- `time` total running time of the algorithm
- `output` path of the output file
- `objective` objective value of the solution

