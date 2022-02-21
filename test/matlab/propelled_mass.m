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

%% Parameters variation
inputStruct.m = 1.0;
inputStruct.dT = 0.1;
inputStruct.T = 2.0;
inputStruct.x0 = 0.1;
inputStruct.cost_multiplier = 1;

experiments = cell(1,12);

initialRelaxed = inputStruct;
initialRelaxed.complementairity = 'Relaxed';
experiments(1:4) = {initialRelaxed};

experiments{1}.epsilon_relaxed = 0.004;
experiments{2}.epsilon_relaxed = 0.008;
experiments{3}.epsilon_relaxed = 0.012;
experiments{4}.epsilon_relaxed = 0.016;

initialDynamical = inputStruct;
initialDynamical.complementairity = 'Dynamical';
experiments(5:8) = {initialDynamical};

experiments{5}.K_dynamic              = 20;
experiments{5}.epsilon_dynamical = 0.05;

experiments{6}.K_dynamic               = 20;
experiments{6}.epsilon_dynamical = 0.1;

experiments{7}.K_dynamic               = 10;
experiments{7}.epsilon_dynamical = 0.05;

experiments{8}.K_dynamic               = 10;
experiments{8}.epsilon_dynamical = 0.1;

initialHyperbolic = inputStruct;
initialHyperbolic.complementairity = 'Hyperbolic';
experiments(9:12) = {initialHyperbolic};

experiments{9}.K_hyperbolic           = 250;
experiments{9}.scaling_hyperbolic = 500;

experiments{10}.K_hyperbolic           = 125;
experiments{10}.scaling_hyperbolic = 500;

experiments{11}.K_hyperbolic           = 250;
experiments{11}.scaling_hyperbolic = 400;

experiments{12}.K_hyperbolic           = 125;
experiments{12}.scaling_hyperbolic = 400;

experimentsResults = struct();

for comp_cell = complementarities
    experimentsResults.(comp_cell{:})  = struct();
    experimentsResults.(comp_cell{:}).elapsedTimes = [];
    experimentsResults.(comp_cell{:}).accuracy = [];
end


for exp = experiments
    [position, velocity, force, propeller, forceDerivative, t, costValue, elapsedTime, freeFalling, expectedForce] = solve_propelled_mass(exp{:});
    experimentsResults.(exp{:}.complementairity).elapsedTimes = [experimentsResults.(exp{:}.complementairity).elapsedTimes elapsedTime];
    experimentsResults.(exp{:}.complementairity).accuracy = [experimentsResults.(exp{:}.complementairity).accuracy mean(force .* position)];
end

figure

for comp_cell = complementarities
    hold on
    scatter(experimentsResults.(comp_cell{:}).elapsedTimes, ...
        experimentsResults.(comp_cell{:}).accuracy, 48, ...
        'filled','s', 'DisplayName', comp_cell{:});
    legend('-DynamicLegend');
end
x_label = xlabel('Average Computational Time (s)');
set(x_label, 'Interpreter', 'latex');
set(x_label, 'FontSize', 16);
y_label = ylabel('Accuracy $p_z \cdot f_z$');
set(y_label,'Interpreter','latex');
set(y_label,'FontSize', 16);

%% initial position take 2
initial_positions = 0.05 : 0.02 : 0.15;
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

