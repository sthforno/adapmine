## Dataset

The full dataset can be download from the original Stanford Large Network Dataset Collection, which can be found at https://snap.stanford.edu/data/.

## Environments

System: Linux

Compilers: gcc@12(minimum 10)

Other dependencies:

+ openmpi (> 4.1)
+ cmake (> 3.21)

## Build

In the root directory of the project:

```bash
mkdir build && cd build
cmake ..
make -j
```

## Usage

in `build/` directory:

`./bin/approx <graph_file> <pattern_file>`



