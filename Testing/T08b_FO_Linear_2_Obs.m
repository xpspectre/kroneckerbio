% Test FO with simple linear model and 2 outputs
%   Model: 0 -> x with 0-th order rate k
%   The analytic solution to this model is x[t] = kt
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
    m = AddSeed(m, 'x1_0', 1);
    m = AddSeed(m, 'x2_0', 2);
    m = AddState(m, 'x1', 'v', 'x1_0');
    m = AddState(m, 'x2', 'v', 'x2_0');
    m = AddParameter(m, 'k1', 0.9);
    m = AddParameter(m, 'k2', 1.0);
    m = AddParameter(m, 'k3', 1.1);
    m = AddReaction(m, 'Reaction1', {}, {'x1'}, 'k1');
    m = AddReaction(m, 'Reaction2', {'x1'}, {'x2'}, '(k2+k3)*x1');
    m = AddOutput(m, 'y1', 'x1');
    m = AddOutput(m, 'y2', 'x2');
    
    Omega = [0.1, 0.01; 0.01, 0.1];
    sigmas = [0.05, 0.01, 0.04];
    Sigma = [sigmas(1), sigmas(2); sigmas(2), sigmas(3)];
    
    m = AddEta(m, 'k1');
    m = AddEta(m, 'k2');
    m = AddOmega(m, {'k1', 'k2'}, Omega);
    
    m = AddErrorModel(m, {'y1','y2'}, {'sigma__y1','';'sigma__y2__y1','sigma__y2'}, {'sigma__y1','sigma__y2__y1','sigma__y2'}, sigmas); % sigma param names are arbitrary
    
    m = FinalizeModel(m);
    
    % Hack this in for now
    thetaNames = {m.Parameters.Name};
    fOmega = calcOmega({'omega__k1__k1','';'omega__k2__k1','omega__k2__k2'}, thetaNames); % use names from AddOmega - TODO: generalize this
    
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
    save('FO_linear_data_2_obs.mat', 'm', 'fOmega', 'times', 'measurements')
    
else
    %% Load data
    load('FO_linear_data_2_obs.mat')
end

%% Test FO method
% Supply analytic solutions and sensitivities
% k = m.k(1);
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
    obs = observationFOCEI({'y1','y2'}, times, 'FOCEI', ['Obs' num2str(i)]);
    obj = obs.Objective(measurements{i}, fOmega);
    
    fit.addFitConditionData(obj, con);
end
opts = [];
opts.Verbose = 1;
opts.ComputeSensPlusOne = true; % NLME require (n+1)-th order sensitivities; Note: calls computeObjGrad, which will try to compute dGdk as well, which isn't needed
opts.Normalized = false; % Simple system has bounded params w/ known ranges
opts.UseAdjoint = false; % Can't use adjoint (for now) since dy/dT sensitivities not calculated
opts.ComplexStep = false; % probably not compatible

G = ObjectiveValue(fit, opts);

[D, G] = ObjectiveGradient(fit, opts);

[Df, Gf] = FiniteObjectiveGradient(fit, opts);

