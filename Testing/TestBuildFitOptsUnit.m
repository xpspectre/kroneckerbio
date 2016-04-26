clear; close all; clc
rng('default');

%% Test set default opts
opts = BuildFitOpts();

%% Test model addition
% Add 1st con
m = equilibrium_model;
con1 = experimentInitialValue(m, [], [], [], 'InitialValueExperiment1');
outputList   = [1    1    2    2    3    3  ]'; % index of the output this measurement refers to
timesList    = [0.1  1    0.1  1    0.1  1  ]';
measurements = [0.6  0.4  1.5  1.3  0.4  0.6]';
sd = sdLinear(0.05, 0.1);
obs = observationLinearWeightedSumOfSquares(outputList, timesList, sd, 'DefaultObservation1');
obj1 = obs.Objective(measurements);

newopts = [];

opts = BuildFitOpts(opts, m, con1, obj1, newopts);

% Add 2nd con with same model as con1 and same UseParams, different seeds
con2 = experimentInitialValue(m, [1.5;1.5;0], [], [], 'InitialValueExperiment2');
outputList   = [1    1    2    2    3    3  ]'; % index of the output this measurement refers to
timesList    = [0.1  1    0.1  1    0.1  1  ]';
measurements = [0.6  0.4  1.5  1.3  0.4  0.6]';
sd = sdLinear(0.05, 0.1);
obs = observationLinearWeightedSumOfSquares(outputList, timesList, sd, 'DefaultObservation2');
obj2 = obs.Objective(measurements);
newopts2 = [];

opts = BuildFitOpts(opts, m, con2, obj2, newopts2);

% Add 3rd con with same model as con1 and different UseParams
con3 = experimentInitialValue(m, [1.6;1.6;0], [], [], 'InitialValueExperiment3');
outputList   = [1    1    2    2    3    3  ]'; % index of the output this measurement refers to
timesList    = [0.1  1    0.1  1    0.1  1  ]';
measurements = [0.6  0.4  1.5  1.3  0.4  0.6]';
sd = sdLinear(0.05, 0.1);
obs = observationLinearWeightedSumOfSquares(outputList, timesList, sd, 'DefaultObservation3');
obj3 = obs.Objective(measurements);
newopts3 = [];
newopts3.UseParams = [1;2]; % 2nd param is independent

opts = BuildFitOpts(opts, m, con3, obj3, newopts3);

% Add 4th con with same model as con1 and different UseParams and 2 obj funs
con4 = experimentInitialValue(m, [1.7;1.7;0], [], [], 'InitialValueExperiment4');
outputList   = [1    1    2    2    3    3  ]'; % index of the output this measurement refers to
timesList    = [0.1  1    0.1  1    0.1  1  ]';
measurements = [0.6  0.4  1.5  1.3  0.4  0.6]';
sd = sdLinear(0.05, 0.1);
obs = observationLinearWeightedSumOfSquares(outputList, timesList, sd, 'DefaultObservation4_1');
obj4_1 = obs.Objective(measurements);
measurements = [0.6  0.4  1.5  1.3  0.4  0.7]';
obs = observationLinearWeightedSumOfSquares(outputList, timesList, sd, 'DefaultObservation4_1');
obj4_2 = obs.Objective(measurements);
newopts4 = [];
newopts4.UseParams = [1;3]; % 2nd param is independent

opts = BuildFitOpts(opts, m, con4, [obj4_1; obj4_2], newopts4);

% Test getting and setting params, bounds, etc
opts.fit.T2Tlocal([4;5;6;7])

% Test running the fit
con = [con1; con2; con3; con4];
obj = [obj1; obj2; obj3; obj4_1; obj4_2];
fit = FitObjective(m, con, obj, opts);


%% Test with multiple model types
% All parameters must have unique indices representing their fit
opts_m = [];
opts_m = BuildFitOpts(opts_m, struct('UseIndependentModels', true));

% Add 1st con
m = struct('Name', 'TestModel1', 'nk', 2);
con1 = struct('Name', 'TestCondition1', 'ns', 3, 'nq', 1, 'nh', 1);
obj1 = struct('Name', 'TestObjective1');
newopts = [];
newopts.UseParams = [1;2];
newopts.UseSeeds = [3;4;5];

opts_m = BuildFitOpts(opts_m, m, con1, obj1, newopts);

% Add 2nd con with different model
m2 = struct('Name', 'TestModel2', 'nk', 3);
con2 = struct('Name', 'TestCondition2', 'ns', 3, 'nq', 1, 'nh', 1);
obj2 = struct('Name', 'TestObjective2');
newopts2 = [];
newopts2.UseParams = [1;2;6];
newopts2.UseSeeds = [7;8;9];

opts_m = BuildFitOpts(opts_m, m2, con2, obj2, newopts2);