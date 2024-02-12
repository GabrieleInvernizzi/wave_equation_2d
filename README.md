# 2D wave equation solver

The project is comprised of two executables, the serial and the parallel implementations. 
The parallel one is built on top of OpenMPI.

## Build
Before trying to build the project ensure to have installed the following dependencies:
- gcc
- make
- OpenMPI (5.0.1)

Next, to download and build the project run the following commands:
```bash
git clone https://github.com/GabrieleInvernizzi/wave_equation_2d.git
cd wave_equation_2d
make
```
This will create two executables in the `./bin` directory:
- `wave_eq_s`: the serial implementation
- `wave_eq_p`: the parallel implementation

There are some make flags that can be setted to tailor compilation, for more info look at the `Makefile`.

## Usage
The serial implementation can be run as any other program.
```bash
./wave_eq_s
```

The parallel implementation must be run as an OpenMPI parallel application, there are multiple ways to do so.
The simplest one is by using `mpirun`.
```bash
mpirun -n 5 ./wave_eq_p
```

To know more about all the options that can be used with both implementations run them with the `help` flag.
```bash
./wave_eq_s --help
```

Next, the simulation file can be encoded as a video using the script `sim2video`.
```bash
pipenv run python3 ./sim2video.py ../out.sim
```

## Example video

https://github.com/GabrieleInvernizzi/wave_equation_2d/assets/42583079/0bcc8ac9-e331-4855-b68e-4ceff9421194

