# Eglin Store Separation Trajectory Prediction

## Introduction
This project involves the development of a Surrogate Model for Eglin Store Separation Trajectory Prediction from Wing using Gaussian Process Regression (GPR). The goal is to predict the forces and moments on a store during separation and simulate its trajectory using a GPR model.

## Code Files Description

- SimData.m
The SimData class is used to estimate forces and moments on the store during separation at each time instant using pressure data on the store at each timestep. It uses a 6DOF file from the simulation to estimate the state variables like effective AOA, sideslip, body angular rates, etc.

- RegressionModel.m
The RegressionModel class uses the simulation data from SimData class instances and builds a GPR model for body force and moment prediction. This model is trained using the processed simulation data and is used to predict the forces and moments on the store.

- sixdofsim.m
The sixdofsim.m script is used for rigid body dynamics (RBD) simulation using the GPR model to predict the store separation trajectory. It loads the GPR model, initializes the simulation parameters, and solves the ODE to simulate the store's trajectory.

- simData.mat
It contains already generated store state and force data using SimData class

## Usage

Ensure all required files and data are in place.
Run the sixdofsim.m script to perform the simulation and predict the store separation trajectory.

## Dependencies

MATLAB
Gaussian Process Regression (GPR) Toolbox


