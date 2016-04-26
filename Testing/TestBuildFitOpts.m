clear; close all; clc
rng('default');

%% Add global opts
opts = [];

% newopts = [];
% newopts.Normalized = false;
% 
% opts = BuildFitOpts(opts, newopts);

%% Add condition to fit
m = equilibrium_model;

con = experimentInitialValue(m, [], [], [], 'InitialValueExperiment');

outputList   = [1    1    2    2    3    3  ]'; % index of the output this measurement refers to
timesList    = [0.1  1    0.1  1    0.1  1  ]';
measurements = [0.6  0.4  1.5  1.3  0.4  0.6]';

% Simple error model with a constant and proportional term
sd = sdLinear(0.05, 0.1);

% Create observation scheme (specifies output and time to observe simulation)
% and objective function (specifies data for comparison)
obs = observationLinearWeightedSumOfSquares(outputList, timesList, sd, 'DefaultObservation');
obj = obs.Objective(measurements);

newopts = [];
newopts.UseParams = [1;1];

opts = BuildFitOpts(opts, m, con, obj, newopts);

%% Add 2nd condition
% opts = BuildFitOpts(opts, m, con, obj, newopts);