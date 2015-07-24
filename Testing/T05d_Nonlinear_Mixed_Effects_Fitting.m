% NONMEM-like mixed effects modeling
clear; close all; clc

%% Initialize Fit Object
fit = FitObject('T05d_NLME_Fit');

%% Load data
theoph = load_theo();

nIds = height(theoph);
% nIds = 3; % subset of patients for now

%% Load model
% Model contains nominal dose of 320 mg
m = theo_model();
% load('theo_model.mat') % preload model for speed after building

fit.addModel(m);

%% Setup conditions for each patient
% Each patient needs a different condition because covariates require each
% patient's ODE system be integrated independently
method = 'FO'; % simplest approx
for i = 1:nIds
    times = theoph(i,:).Time{1};
    concs = theoph(i,:).Conc{1};
%     dose = theoph(i,:).Dose; % real covariate
    dose = 320; % fix all patients for testing
    
    % Specify condition
    con = experimentInitialValue(m, [dose; 0], [], [], ['TheoCon' num2str(i)]); % modified ICs
    
    % Specify objective function
    obs = observationNLME('Cp', times, method, ['TheoObs' num2str(i)]);
    obj = obs.Objective(concs);
    
    % Specify table of parameter, condition, and fitting options
    Name = {'ka', 'kel', 'V'}';
    Type = {'k',  'k',   'k'}';
    LB   = [0,    0,     0]';
    UB   = [3     0.3    50]';
    Use  = [1     1      1]';
    
    opts = [];
    opts.ParamSpec = table(Name, Type, LB, UB, Use);
    
    % Add patient to fit
    fit.addFitConditionData(obj, con, opts);
end

%% Plot data
% sim = SimulateSystem(m, con, 25); % simulate base system
% 
% figure
% hold on
% % Model
% t = linspace(0,25,100)';
% plot(t, sim.y(t,2));
% % Data
% for i = 1:nIds
%     plot(theoph(i,:).Time{1}, theoph(i,:).Conc{1}, '+')
% end
% xlabel('Time (hr)')
% ylabel('Cp (mg/L)')
% legend('Model', 'Data')
% hold off

%% Setup fit
opts.Verbose = 2;
opts.Normalized = false; % Simple system has bounded params w/ known ranges
opts.UseAdjoint = false; % can't use adjoint (for now) since dy/dT sensitivities not calculated

%% Test objective value
% G = ObjectiveValue(fit, opts); % doesn't work because NLME requires gradients even for G calculation
[~, G] = ObjectiveGradient(fit, opts); % will get 33 warnings: 3 patients, 11 timepoints each

%% Fit
fitOut = FitObjective(fit, opts);
