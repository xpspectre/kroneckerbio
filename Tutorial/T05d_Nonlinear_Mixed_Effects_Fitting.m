%% T05d_ Nonlinear Mixed Effects Fitting
% NONMEM-like mixed effects modeling

%% Load data
addpath([fileparts(mfilename('fullpath')) '/../Testing']);
theoph = load_theo();

% Number of patients to fit.
% nIds = height(theoph); % all patients
nIds = 3; % subset of patients for faster testing

% Visualize dataset
figure
hold on
for i = 1:nIds
    plot(theoph.Time{i}, theoph.Conc{i})
end
hold off
xlabel('Time')
ylabel('Conc')
title(sprintf('Theophylline dataset for %g Patients', nIds))

%% Load model
buildModel = true; % Model building can take a while so save and load subsequent tries

if buildModel
    m = theo_model;
    m_nlme = theo_model('nlme');
    save('theo_model_.mat', 'm', 'm_nlme');
else
    load('theo_model_.mat');
end

%% Run basic fit of averaged data
% To get sane starting values of parameter population avg values
% Get time range
tmax = 0;
for i = 1:nIds
    tmax_i = theoph(i,:).Time{1}(end);
    if tmax_i > tmax
        tmax = tmax_i;
    end
end
nTimes = 20;
times = linspace(0, tmax, nTimes)';

% Get average concentration and dose data by interpolation
concs = zeros(nTimes,nIds);
doses = zeros(1,nIds);
for i = 1:nIds
    times_i = theoph(i,:).Time{1};
    concs_i = theoph(i,:).Conc{1};
    concs(:,i) = interp1(times_i, concs_i, times, 'pchip', 'extrap'); % small amount of extrapolation won't go negative
    doses(i) = theoph(i,:).Dose;
end
concAvg = mean(concs, 2);
concStd = std(concs, [], 2);
dose = mean(doses);

sd = stdError(times, concStd); % Fix variance to measured variance

con = experimentInitialValue(m, [dose; 0], [], [], 'TheoConAvg'); % modified ICs

CpInd = find(ismember({m.Outputs.Name}, 'Cp'));
outputList = repmat(CpInd, nTimes, 1);

obs = observationLinearWeightedSumOfSquares(outputList, times, sd, 'TheoObsAvg');
obj = obs.Objective(concAvg);

opts = [];
opts.TolOptim = 1e-3; % large tolerance to run fast for test

mFit = FitObjective(m, con, obj);

% Display avg fit
sim = SimulateSystem(mFit, con, tmax);
t = linspace(0, tmax, 100)';
theo_avg_fit = sim.y(t, CpInd)';

figure
hold on
% Model
plot(t, theo_avg_fit);
% Data
for i = 1:nIds
    plot(theoph(i,:).Time{1}, theoph(i,:).Conc{1}, '+')
end
xlabel('Time (hr)')
ylabel('Cp (mg/L)')
legend('Model', 'Data')
title('Comparison of Avg Fit to Data')
hold off

%% Setup conditions for each patient
% Each patient needs a different condition because covariates require each
% patient's ODE system be integrated independently
fitNlme = FitObject('Theo_NLME_Fit');
fitNlme.addModel(m_nlme);

for i = 1:nIds
    times = theoph(i,:).Time{1}';
    concs = theoph(i,:).Conc{1}';
    
    dose = theoph(i,:).Dose; % real covariate
%     dose = 320; % fix all patients for testing
    
    % Specify condition
    con = experimentInitialValue(m, [dose; 0], [], [], ['TheoCon' num2str(i)]); % modified ICs
    
    % Specify objective function
    obs = observationFocei('Cp', times, 'FO', ['TheoObs' num2str(i)]);
    obj = obs.Objective(concs);
    
    % Specify table of parameter, condition, and fitting options
    opts = [];
    
    % Optionally specify params to fit and bounds. By default, only params
    %   specified here with Use not equal to 0 are fit. Commented out to fit
    %   everything.
%     Name = {'ka', 'kel', 'V'}';
%     Type = {'k',  'k',   'k'}';
%     LB   = [ 0,    0,     0 ]';
%     UB   = [ 5,    0.3,   50]';
%     Use  = [ 1,    1,     1 ]';
%     opts.ParamSpec = table(Name, Type, LB, UB, Use, 'RowNames', Name);
    
    % Add patient to fit
    fitNlme.addFitConditionData(obj, con, opts);
end

%% Setup fit
opts.Verbose = 2;
opts.ComputeSensPlusOne = true; % NLME require (n+1)-th order sensitivities; Note: calls computeObjGrad, which will try to compute dGdk as well, which isn't needed
opts.Normalized = false; % Simple system has bounded params w/ known ranges
opts.UseAdjoint = false; % Can't use adjoint (for now) since dy/dT sensitivities not calculated
opts.ComplexStep = false; % probably not compatible
opts.TolOptim = 1e-1; % large tolerance to run fast for test

%% Test objective value
G = ObjectiveValue(fitNlme, opts);

%% Fit
fitNlmeOut = FitObjective(fitNlme, opts);

%% Display average fit result
% TODO: Display individual fits, calculating final etas
con = fitNlmeOut.Conditions(i);
mFitNlme = fitNlmeOut.Models(con.ParentModelIdx);
sim = SimulateSystem(mFitNlme, con, tmax);
theo_nlme_fit = sim.y(t,CpInd)';

figure
hold on
% Model
plot(t, theo_nlme_fit);
% Data
for i = 1:nIds
    plot(theoph(i,:).Time{1}, theoph(i,:).Conc{1}, '+')
end
xlabel('Time (hr)')
ylabel('Cp (mg/L)')
legend('Avg Model', 'Data')
title('Comparison of NLME Fit to Data')
hold off
