% Run a simple simulation and plot the results.
% Do some additional analysis: sensitivity and curvature calculation, and linear
%   noise approximation to make sure they work.

%% Load equilibrium experiment A + B <-> C
m = LoadModelMassAction('Equilibrium.txt');
% m = equilibrium_model_analytic;

%% Construct experiment
con = experimentInitialValue(m, [], [], [], 'InitialValueExperiment');
con = ExperimentInitialValue(m, [], [], [], 'InitialValueExperiment');

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

%% Modify stuff and test
tic
for i = 1:1000
    m = m.Update(rand(2,1));
    sim1 = SimulateSystem(m, con, tF);
end
toc

m = m.Update([1,0.5]');

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