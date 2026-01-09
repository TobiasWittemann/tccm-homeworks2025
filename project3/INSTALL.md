## Usage
### Compilation
```bash
gfortran src/MD.f90 -o md
```
### Simulation
Execute the binary and follow the interactive prompts to enter the md configurations:
```bash
./md
```
User input:
- Input filename: Path to your initial coordinates (Press **Enter** for default: inp.txt).

- Number of steps: Total iterations (Press **Enter** for default: 1000).

- Time step (dt): Integration step in ps (Press **Enter** for default: 0.02).

### Output
The simulation generates a `traj.xyz` file which is recorded every 10 steps.
