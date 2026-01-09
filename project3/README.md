# MD: Argon Molecular Dynamics Simulation
A Molecular Dynamics (MD) simulation program written in Fortran. This program simulates a system of Argon atoms interacting via the **Lennard-Jones potential**, utilizing the **Velocity Verlet** algorithm for numerical integration of the equations of motion.

Run the python script (track_energy.py) to track the E(tot) in function of time steps it in the same directory where traj.xyz will be located or provide a path to traj.xyz when running.

# Directory structure
```text
project3/
├── src/
│   ├── MD.f90
│   └── track_energy.py
├── tests/
│   ├── inp.txt   # Test input
│   └── traj.xyz  # Expected output when using inp.txt
├── AUTHORS
├── INSTALL.md
├── LICENSE
├── README.md
└── dynamics.pdf
