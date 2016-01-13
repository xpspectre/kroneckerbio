% Test FO with simple linear model
%   Model: 0 -> x with 0-th order rate k
%   The analytic solution to this model is x[t] = kt
%       including the initial conditions is x[t] = x0 + kt
%   Inter-individual variability is for k, with the form k = k0*exp(eta_k)
%   Intra-individual variability is for x, with the form eps ~ N(0, sigma_k^2)

clear; close all; clc
rng('default');

% Whether to generate new data for load model and simulated data
simNew = false;

if simNew
    %% Build model
    m = InitializeModelAnalytic('LinearModel');
    m = AddCompartment(m, 'v', 3, 1);
    m = AddSeed(m, 'x0', 1);
    m = AddState(m, 'x', 'v', 'x0');
    m = AddParameter(m, 'k', 1.5);
    m = AddReaction(m, 'Production', {}, {'x'}, 'k');
    m = AddOutput(m, 'y', 'x');
    
    Omega = 0.2; % Omega specified as omegas^2 (for now)
    Sigma = sqrt(0.05); % Sigma specified as parameters in R
    
    m = AddEta(m, 'k');
    m = AddOmega(m, {'k'}, Omega);
    
    m = AddErrorModel(m, {'y'}, {'sigma__y^2'}, {'sigma__y'}, Sigma);
    
    m = FinalizeModel(m);
    
    %% Generate simulated data
    n = 3; % number of patients
    nPoints = 5; % number of measurements/patient
    tf = 10; % final time
    times = linspace(0, tf, nPoints)';
    
    measurements = zeros(nPoints, n);
    for i = 1:n
        % Simulate 1 patient
        % Draw from variability distributions
        eta_k = normrnd(0, sqrt(Omega));
%         eta_k = 0; % just replace with 0 for a test case
        err = normrnd(0, sqrt(Sigma), nPoints, 1); % generate
%         err = repmat(0.1, nPoints, 1); % just replace with 0.1 for a test case
        
        % Make model for patient
        mi = m;
        k = mi.k;
        k(1) = k(1)*exp(eta_k); % hardcode in inter-individual variability for this case
        mi = mi.Update(k); % remember to reassign since this is a struct/by value
        
        con = experimentInitialValue(mi, [], [], [], ['Con' num2str(i)]);
        
        % Simulate patient
        sim = SimulateSystem(mi, con, tf);
        
        % Get measurements
        measurement = sim.x(times,1)' + err; % hardcode in 1st state = output + intra-indivual variability
        measurements(:,i) = measurement;
        
    end
    
    %% Save data
    save('FO_linear_data.mat', 'm', 'times', 'measurements')
    
else
    %% Load data
    load('FO_linear_data.mat')
end

%% Test FO method
% Supply analytic solutions and sensitivities
n = size(measurements, 2);
tf = times(end);

% Do basic simulation with nominal params (0 inter- and intra-individual variability)
con = experimentInitialValue(m, [], [], [], 'Con');
sim = SimulateSystem(m, con, tf);
exact = sim.x(times, 1)';

% Make fitting object
fit = FitObject('Fit_FO_Linear');
fit.addModel(m);
for i = 1:n
    con = experimentInitialValue(m, [], [], [], ['Con' num2str(i)]);
    
    obs = observationFocei('y', times, 'FO', ['Obs' num2str(i)]);
    obj = obs.Objective(measurements(:,i)'); % measurements must be specified as a nOutput x ntimes matrix
    
    fit.addFitConditionData(obj, con);
end
opts = [];
opts.Verbose = 2; % 2 to show inner optimization steps
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
% Specify parameters to fit; default is to fit all params
% TODO: Right now, all params are fit; implement specifying which params
%   are fit, and automatic specification of corresponding covariances
% TODO: In FOCEI, not fitting an eta means setting it to 0 always
% Default theta bounds are [0, Inf)
% Default eta bounds are (-Inf, Inf), but aren't explicitly fit
% Note: Right now, params that make up covariance matrices aren't enforced
% to make valid covariance matrices
%   But these would probably cause the model to be ill-conditioned - how
%   does this appear in the optimization?
opts.Approximation = 'FO';
% opts.Approximation = 'FOCEI';
opts.MaxStepSize             = 0.1; % smaller (jumps around too much with normal FitObjective value of 1)
opts.MaxFunEvals             = 20000; % more fun evals/step 
opts.MaxIter                 = 1000;

innerOpts = [];
opts.InnerOpts = innerOpts;

[fitOut, G, D] = FitObjectiveNlme(fit, opts);

% FO final params values for default single participant:
% k = [1.51668703306234]
% omega__k__k = [0.0598858830220903]
% sigma__y = [0.165020337068679]

% FO final params values for default 3 participants:
% k = [1.92645269775046]
% omega__k__k = [0.00411617548244264]
% sigma__y = [0.373219359558438]

