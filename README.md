# README for SchrodingerSolver Repository

## Overview
The `SchrodingerSolver` repository contains a MATLAB function designed to solve the Schrödinger equation for different potential wells commonly encountered in quantum mechanics. This MATLAB script is aimed at students and researchers in physics, particularly those working in the field of quantum mechanics.

## Features
- **Flexible Potential Well Selection**: Users can choose from three different potential wells:
  - Infinite Square Well
  - Finite Square Well
  - Harmonic Oscillator
- **User Input for Parameters**: The script allows users to input parameters like well width, effective mass, and the number of grid points for calculations.
- **Visualization of Results**: It generates plots for eigenenergies and eigenfunctions, providing a visual understanding of quantum states within the selected potential well.
- **Theoretical and Calculated Energy Comparison**: For certain potential wells, the script calculates and displays both theoretical and computed energies, enabling an analysis of the solution's accuracy.

## Data Handling
The data is handled using standard physics constants and user-provided inputs. The script dynamically adjusts calculations based on the selected potential type and user inputs.

## Usage
Run the `schrodingerSolver` function in MATLAB. Follow the prompts to input the desired parameters for the potential well and the particle. The script will compute and display the eigenenergies and eigenfunctions for the given setup.

## Technical Description
The script uses fundamental constants like reduced Planck's constant, elementary charge, and electron mass for its calculations. It employs numerical methods to solve the eigenvalue problem of the Hamiltonian matrix, constructed using finite differences. The potential V is set according to the user's choice of well type, and eigenenergies are calculated in electron volts (eV). The script also includes a helper function `getWellType` to convert the well type numeric input into a readable string format.

## Note
This project is ideal for educational and research purposes and offers a practical approach to understanding quantum mechanics solutions.

