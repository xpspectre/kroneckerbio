% Run a simple simulation and plot the results.
% Do some additional analysis: sensitivity and curvature calculation, and linear
%   noise approximation to make sure they work.

%% Load equilibrium experiment A + B <-> C
m = LoadModelMassAction('Equilibrium.txt');

%% Construct experiment
con = experimentInitialValue(m, [], [], [], 'InitialValueExperiment');

%% Simulate
tF = 1; % final time
sim1 = SimulateSystem(m, con, tF);

% Plot result
figure
times = linspace(0, tF, 100);
plot(times, sim1.x(times))
legend('A','B','C')
xlabel('Time')
ylabel('Amount')

%% Simulate Sensitivities
sim2 = SimulateSensitivity(m, con, tF);

%% Simulate Curvature
opts.UseSeeds = [];
sim3 = SimulateCurvature(m, con, tF, opts);

%% Simulate Linear Noise Approximation
obs = observationAll(tF);
sim4 = SimulateLna(m, con, obs);