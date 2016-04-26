%% T05c Fit Multiple Models
% Test/tutorial script for multiple distinct model/topology fits

%% Initialize fit object
opts = BuildFitOpts;
opts.Verbose = 2;
opts.TolOptim = 0.01;
opts.Normalized = false;
opts.UseIndependentModels = true; % allows multiple model types and relationships between fit params

%% Construct equilibrium experiment A + B <-> C with seeds
addpath([fileparts(mfilename('fullpath')) '/../Testing']);
m1 = equilibrium_model;

%% Construct 2nd model that requires dimerization of A
m2 = equilibrium_dimer_model;

%% Set up experimental conditions and generate some test data for each model
tF = 1;
nTimes = 5;
times = linspace(0, tF, nTimes)';
outputs = {'A','B','C'};
sd = sdLinear(0.05, 0.1);
    
% Sample model 1 experiment and data
con1 = experimentInitialValue(m1, [1,2,0], [], [], 'Condition1');

[outputsList, timesList, measurementsList] = generateTestData(m1, con1, times, outputs, sd);
measurements1 = reshape(measurementsList, nTimes, 3);

obs = observationLinearWeightedSumOfSquares(outputsList, timesList, sd, 'Observation1');
obj1 = obs.Objective(measurementsList);

newopts = [];
newopts.UseParams = [1,2]';
newopts.UseSeeds = [1,2,0]';
newopts.ParamLowerBound = [1e-2 1e-2];
newopts.ParamUpperBound = [1e2  1e2];

opts = BuildFitOpts(opts, m1, con1, obj1, newopts);

% Sample model 2 experiment and data
con2 = experimentInitialValue(m2, [2,2,0,0], [], [], 'Condition2');

[outputsList, timesList, measurementsList] = generateTestData(m2, con2, times, outputs, sd);
measurements2 = reshape(measurementsList, nTimes, 3);

obs = observationLinearWeightedSumOfSquares(outputsList, timesList, sd, 'Observation2');
obj2 = obs.Objective(measurementsList);

newopts = [];
newopts.UseParams = [3,2,4,5]';
newopts.UseSeeds = [3,4,0,0]';
newopts.ParamLowerBound = [1e-2 1e-2 1e-2 1e-2];
newopts.ParamUpperBound = [1e2  1e2  1e2  1e2];

opts = BuildFitOpts(opts, m2, con2, obj2, newopts);

%% Run fit
m = [m1;m2];
con = [con1;con2];
obj = [obj1;obj2];
[mfit, confit, G] = FitObjective(m, con, obj, opts);

%% Display fit results
timesFine = linspace(0, tF, 100)';
simFit1 = SimulateSystem(mfit(1), confit(1), tF);
simFit2 = SimulateSystem(mfit(2), confit(2), tF);

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


    