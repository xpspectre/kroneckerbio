% Test/tutorial script for multiple distinct model/topology fits
clear; close all; clc
rng('default')

%% Initialize fit object
opts = [];
opts.Verbose = 2;
opts.TolOptim = 0.01;
opts.Normalized = false;
opts.ModelsShareParams = true; % allows multiple model types and relationships between fit params
fit = FitObject('T09c_Fit_Multiple_Models', opts);

%% Construct equilibrium experiment A + B <-> C with seeds
m1 = equilibrium_model;
fit.addModel(m1);

%% Construct 2nd model that requires dimerization of A
m2 = equilibrium_dimer_model;
fit.addModel(m2);

%% Set up experimental conditions and generate some test data for each model
tF = 1;
nTimes = 5;
times = linspace(0, tF, nTimes)';
outputs = {'A','B','C'};
sd = sdLinear(0.05, 0.1);
    
% Sample model 1 experiment and data
con = experimentInitialValue(m1, [1,2,0], [], [], 'Condition1');

[outputsList, timesList, measurementsList] = generateTestData(m1, con, times, outputs, sd);
measurements1 = reshape(measurementsList, nTimes, 3);

obs = observationLinearWeightedSumOfSquares(outputsList, timesList, sd, 'Observation1');
obj = obs.Objective(measurementsList);

opts = [];
opts.UseParams = [1,2]';
opts.UseSeeds = [1,2,0]';
opts.ParamLowerBound = [1e-2 1e-2];
opts.ParamUpperBound = [1e2  1e2];

fit.addFitConditionData(obj, con, opts);

% Sample model 2 experiment and data
con = experimentInitialValue(m2, [2,2,0,0], [], [], 'Condition2');

[outputsList, timesList, measurementsList] = generateTestData(m2, con, times, outputs, sd);
measurements2 = reshape(measurementsList, nTimes, 3);

obs = observationLinearWeightedSumOfSquares(outputsList, timesList, sd, 'Observation2');
obj = obs.Objective(measurementsList);

opts = [];
opts.UseParams = [3,2,4,5]';
opts.UseSeeds = [3,4,0,0]';
opts.ParamLowerBound = [1e-2 1e-2 1e-2 1e-2];
opts.ParamUpperBound = [1e2  1e2  1e2  1e2];

fit.addFitConditionData(obj, con, opts);

%% Run fit
fitOut = FitObjective(fit);

%% Display fit results
timesFine = linspace(0, tF, 100)';
simFit1 = SimulateSystem(fitOut.Models(1), fitOut.Conditions(1), tF);
simFit2 = SimulateSystem(fitOut.Models(2), fitOut.Conditions(2), tF);

figure
subplot(1,2,1)
hold on
plot(timesFine, simFit1.y(timesFine)')
ax = gca;
ax.ColorOrderIndex = 1;
plot(times, measurements1, '+')
xlabel('Time')
ylabel('Amount')
legend('A','B','C')
title('Model 1')
hold off

subplot(1,2,2)
hold on
plot(timesFine, simFit2.y(timesFine)')
ax = gca;
ax.ColorOrderIndex = 1;
plot(times, measurements2, '+')
xlabel('Time')
legend('A','B','C','A:A')
title('Model 2')
hold off


    