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

labels = {};
for comp_cell = complementarities
    labels = [labels, comp_cell];
end

%% First run
for comp_cell = complementarities
    inputStruct.complementairity = comp_cell{:};

    [position, velocity, force, propeller, forceDerivative, t, costValue, elapsedTime, freeFalling, expectedForce] = solve_propelled_mass(inputStruct);
    plot_propelled_mass_output(inputStruct, position, velocity, force, propeller, forceDerivative, t, costValue, elapsedTime, freeFalling, expectedForce);

end

%% masses
masses = linspace(0.05, 1.05, 20);
masses_result = struct();

for comp_cell = complementarities
    masses_result.(comp_cell{:}) = struct();
    masses_result.(comp_cell{:}).elapsed_times = zeros(1,length(masses));
    masses_result.(comp_cell{:}).complementarity_average = zeros(1,length(masses));
    masses_result.(comp_cell{:}).costValue = zeros(1,length(masses));
end


for i = 1 : length(masses)
    i
    inputStruct.m = masses(i);
    for comp_cell = complementarities
        inputStruct.complementairity = comp_cell{:};

        [position, velocity, force, propeller, forceDerivative, t, costValue, elapsedTime, freeFalling, expectedForce] = solve_propelled_mass(inputStruct);
        masses_result.(comp_cell{:}).elapsed_times(i) = elapsedTime;
        masses_result.(comp_cell{:}).complementarity_average(i) = mean(position .* force);
        masses_result.(comp_cell{:}).costValue(i) = costValue;
    end
end

inputStruct.m = 0.2;

figure
time_elapsed_matrix = [];
for comp_cell = complementarities
    hold on
    plot(masses, masses_result.(comp_cell{:}).elapsed_times)
    time_elapsed_matrix = [time_elapsed_matrix; masses_result.(comp_cell{:}).elapsed_times];
end

hold on
plot(masses, min(time_elapsed_matrix, [], 1), '*');
legend([labels, 'Best Time'])



figure
for comp_cell = complementarities
    hold on
    plot(masses, masses_result.(comp_cell{:}).complementarity_average)
end
legend(labels)

figure
for comp_cell = complementarities
    hold on
    plot(masses, masses_result.(comp_cell{:}).costValue)
end
legend(labels)

%% initial position
initial_positions = linspace(0.05, 0.25, 20);
positions_result = struct();

for comp_cell = complementarities
    positions_result.(comp_cell{:}) = struct();
    positions_result.(comp_cell{:}).elapsed_times = zeros(1,length(initial_positions));
    positions_result.(comp_cell{:}).complementarity_average = zeros(1,length(initial_positions));
    positions_result.(comp_cell{:}).costValue = zeros(1,length(initial_positions));
end


for i = 1 : length(initial_positions)
    i
    inputStruct.x0 = initial_positions(i);
    for comp_cell = complementarities
        inputStruct.complementairity = comp_cell{:};

        [position, velocity, force, propeller, forceDerivative, t, costValue, elapsedTime, freeFalling, expectedForce] = solve_propelled_mass(inputStruct);
        positions_result.(comp_cell{:}).elapsed_times(i) = elapsedTime;
        positions_result.(comp_cell{:}).complementarity_average(i) = mean(position .* force);
        positions_result.(comp_cell{:}).costValue(i) = costValue;
    end
end
inputStruct.x0 = 0.1;

figure
time_elapsed_matrix = [];
for comp_cell = complementarities
    hold on
    plot(initial_positions, positions_result.(comp_cell{:}).elapsed_times)
    time_elapsed_matrix = [time_elapsed_matrix; positions_result.(comp_cell{:}).elapsed_times];
end
hold on
plot(initial_positions, min(time_elapsed_matrix, [], 1), '*');
legend([labels, 'Best Time'])

figure
for comp_cell = complementarities
    hold on
    plot(initial_positions, positions_result.(comp_cell{:}).complementarity_average)
end
legend(labels)

figure
for comp_cell = complementarities
    hold on
    plot(initial_positions, positions_result.(comp_cell{:}).costValue)
end
legend(labels)


