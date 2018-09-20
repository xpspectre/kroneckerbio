%% T08a NLME FO Simple Linear Model 2 Outputs Test
% Test FO with simple linear model and 2 outputs
%   Model: x1 -> x2 with rate k1
%   The analytic/symbolic solution to this model is tracatable:
%       y1 = x1_0*exp(-k1*t);
%       y2 = x1_0 - x1_0*exp(-k1*t) + x2_0;
%   Inter-individual variability is for k1, with the form k1 = k1_0*exp(eta_k1)
%       An additional parameter k2 is added, which should give a bunch of 0's
%   Intra-individual variability is for x1 and x2, with variance and covariance

clear; close all; clc
rng('default');

% Whether to generate new data for load model and simulated data
simNew = true;

if simNew
    %% Build model
    m = InitializeModelAnalytic('LinearModel');
    m = AddCompartment(m, 'v', 3, 1);
    m = AddSeed(m, 'x1_0', 10);
    m = AddSeed(m, 'x2_0', 0.1);
    m = AddState(m, 'x1', 'v', 'x1_0');
    m = AddState(m, 'x2', 'v', 'x2_0');
    m = AddParameter(m, 'k1', 0.9);
    m = AddParameter(m, 'k2', 1.0);
    m = AddReaction(m, 'Reaction1', {'x1'}, {'x2'}, 'k1*x1');
    m = AddOutput(m, 'y1', 'x1');
    m = AddOutput(m, 'y2', 'x2');
    
    Omega = [0.1, 0.01; 0.01, 0.1];
%     sigmas = [0.05, 0.01, 0.04];
    sigmas = [0.05, 0.04];
%     Sigma = [sigmas(1), sigmas(2); sigmas(2), sigmas(3)];
    Sigma = [sigmas(1), 0; 0, sigmas(2)];
    
    m = AddEta(m, 'k1');
    m = AddEta(m, 'k2');
    m = AddOmega(m, {'k1', 'k2'}, Omega);
    
%     m = AddErrorModel(m, {'y1','y2'}, {'sigma__y1^2','';'sigma__y2__y1^2','sigma__y2^2'}, {'sigma__y1','sigma__y2__y1','sigma__y2'}, sigmas); % sigma param names are arbitrary
    m = AddErrorModel(m, {'y1','y2'}, {'sigma__y1^2','';'','sigma__y2^2'}, {'sigma__y1','sigma__y2'}, sigmas); % sigma param names are arbitrary

    m = FinalizeModel(m);
    
    %% Generate simulated data
    n = 1; % number of patients
    nPoints = 5; % number of measurements/patient (for now, both measurements are at the same times)
    tf = 10; % final time
    times = linspace(0, tf, nPoints)';
    
    measurements = cell(n,1);
    for i = 1:n
        % Simulate 1 patient
        % Draw from variability distributions
        Eta = mvnrnd([0, 0], Omega);
        Eps = mvnrnd([0, 0], sqrt(Sigma), nPoints);
        % Simple error terms for testing
%         Eta = [0 0];
%         Eps = repmat(0.1, nPoints, 2);
        
        % Make model for patient
        mi = m;
        k = mi.k;
        k(1) = k(1)*exp(Eta(1)); % hardcode in inter-individual variability for this case
        k(2) = k(2)*exp(Eta(2));
        mi = mi.Update(k); % remember to reassign since this is a struct/by value
        
        con = experimentInitialValue(mi, [], [], [], ['Con' num2str(i)]);
        
        % Simulate patient
        sim = SimulateSystem(mi, con, tf);
        
        % Get measurements
        measurement = sim.x(times,1:2)' + Eps; % hardcode in n-th state = output + intra-indivual variability
        measurements{n} = measurement;
        
    end
    
    %% Save data
    save('FO_linear_data_2_obs.mat', 'm', 'times', 'measurements')
    
else
    %% Load data
    load('FO_linear_data_2_obs.mat')
end

%% Test FO method
% Supply analytic solutions and sensitivities
n = numel(measurements);
tf = times(end);

% Do basic simulation with nominal params (0 inter- and intra-individual variability)
con = experimentInitialValue(m, [], [], [], 'Con');
sim = SimulateSystem(m, con, tf);
exact = sim.x(times, 1)';

% Test G and dG for individual
fit = FitObject('Fit_FO_Linear_2_Obs');
fit.addModel(m);
for i = 1:n
    con = experimentInitialValue(m, [], [], [], ['Con' num2str(i)]);
    
    obs = observationFocei({'y1','y2'}, times, 'FOCEI', ['Obs' num2str(i)]);
    obj = obs.Objective(measurements{i}'); % measurements must be specified as a nOutput x ntimes matrix
    
    fit.addFitConditionData(obj, con);
end
opts = [];
opts.Verbose = 3;
opts.ComputeSensPlusOne = true; % NLME require (n+1)-th order sensitivities; Note: calls computeObjGrad, which will try to compute dGdk as well, which isn't needed
opts.Normalized = false; % Simple system has bounded params w/ known ranges
opts.UseAdjoint = false; % Can't use adjoint (for now) since dy/dT sensitivities not calculated
opts.ComplexStep = false; % probably not compatible

% G = ObjectiveValue(fit, opts)
% 
% [D, G] = ObjectiveGradient(fit, opts)
% 
% [Df, Gf] = FiniteObjectiveGradient(fit, opts)

%% Run fit
% opts.Approximation = 'FO';
opts.Approximation = 'FOCEI';
opts.MaxStepSize             = 0.1; % smaller (jumps around too much with normal FitObjective value of 1)
opts.MaxFunEvals             = 20000; % more fun evals/step 
opts.MaxIter                 = 1000;

innerOpts = [];
opts.InnerOpts = innerOpts;

initialVals = [fit.Models.Parameters.Value]';
[fitOut, G, D] = FitObjectiveNlme(fit, opts);
finalVals = [fitOut.Models.Parameters.Value]'; % fit's values are modified throughout (since it's a handle object) so be careful

% Show final params
fprintf('Initial param values are %g\n', initialVals)
fprintf('Final   param values are %g\n', finalVals)