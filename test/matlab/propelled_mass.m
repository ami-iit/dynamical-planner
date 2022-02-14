close all 
clear all

%% Constants
inputStruct = struct();
inputStruct.m = 0.2;
inputStruct.g = -9.81;

inputStruct.dT = 0.01;
inputStruct.T = 1;

inputStruct.x0 = 0.1;
inputStruct.v0 = 0;
inputStruct.f0 = 0;

inputStruct.cost_multiplier = 100;

inputStruct.M_fdot = 100;
inputStruct.epsilon_relaxed = 0.004;
inputStruct.K_dynamic = 20;
inputStruct.epsilon_dynamical = 0.05;
inputStruct.K_hyperbolic = 250;
inputStruct.scaling_hyperbolic = 500;

complementarities = {'Relaxed', 'Dynamical', 'Hyperbolic'};

for comp_cell = complementarities
    inputStruct.complementairity = comp_cell{:};

    [position, velocity, force, propeller, forceDerivative, t, costValue, elapsedTime, freeFalling, expectedForce] = solve_propelled_mass(inputStruct);
    plot_propelled_mass_output(inputStruct, position, velocity, force, propeller, forceDerivative, t, costValue, elapsedTime, freeFalling, expectedForce);

end


