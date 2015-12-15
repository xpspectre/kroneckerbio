% Test FO with simple linear model
%   Model: 0 -> x with 0-th order rate k
%   The analytic solution to this model is x[t] = kt
%       including the initial conditions is x[t] = x0 + kt
%   Inter-individual variability is for k, with the form k = k0*exp(eta_k)
%   Intra-individual variability is for x, with the form eps ~ N(0, sigma_k^2)

clear; close all; clc
rng('default');

% Whether to generate new data for load model and simulated data
simNew = true;

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
    
    % Hack this in for now
    thetaNames = {m.Parameters.Name};
    fOmega = calcOmega({'omega__k__k'}, thetaNames); % use names from AddOmega - TODO: generalize this
    
    m = FinalizeModel(m);
    
    %% Generate simulated data
    n = 1; % number of patients
    nPoints = 5; % number of measurements/patient
    tf = 10; % final time
    times = linspace(0, tf, nPoints)';
    
    measurements = zeros(nPoints, n);
    for i = 1:n
        % Simulate 1 patient
        % Draw from variability distributions
%         eta_k = normrnd(0, sqrt(Omega));
        eta_k = 0; % just replace with 0 for a test case
%         err = normrnd(0, sqrt(Sigma), nPoints, 1); % generate
        err = repmat(0.1, nPoints, 1); % just replace with 0.1 for a test case
        
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
    save('FO_linear_data.mat', 'm', 'fOmega', 'times', 'measurements')
    
else
    %% Load data
    load('FO_linear_data.mat')
end

%% Test FO method
% Supply analytic solutions and sensitivities
k = m.k(1);
n = size(measurements, 2);
tf = times(end);

% Do basic simulation with nominal params (0 inter- and intra-individual variability)
con = experimentInitialValue(m, [], [], [], 'Con');
sim = SimulateSystem(m, con, tf);
exact = sim.x(times, 1)';

% Test G and dG for individual
fit = FitObject('Fit_FO_Linear');
fit.addModel(m);
for i = 1:n
    obs = observationFOCEI('y', times, 'FOCEI', ['Obs' num2str(i)]);
    obj = obs.Objective(measurements(:,i), fOmega);
    
    fit.addFitConditionData(obj, con);
end
opts = [];
opts.Verbose = 1;
opts.ComputeSensPlusOne = true; % NLME require (n+1)-th order sensitivities; Note: calls computeObjGrad, which will try to compute dGdk as well, which isn't needed
opts.Normalized = false; % Simple system has bounded params w/ known ranges
opts.UseAdjoint = false; % Can't use adjoint (for now) since dy/dT sensitivities not calculated
opts.ComplexStep = false; % probably not compatible

G = ObjectiveValue(fit, opts)

[D, G] = ObjectiveGradient(fit, opts)

[Df, Gf] = FiniteObjectiveGradient(fit, opts)


