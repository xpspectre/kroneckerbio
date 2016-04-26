%% T05b Multiple Participants Fitting
% Test/tutorial script for multiple participants/experiments with
%   datasets/objective functions, where some parameters are shared between
%   participants and others are independent.

%% Initialize fit options
opts = BuildFitOpts;

%% Construct equilibrium experiment A + B <-> C with seeds
addpath([fileparts(mfilename('fullpath')) '/../Testing']);
m = equilibrium_model;

%% Set up experimental conditions and generate some test data for multi-"participant" fitting
tF = 1;
times = linspace(0, tF, 5)';
outputs = {'A','B','C'};
sd = sdLinear(0.05, 0.1);

nCon = 3;
con = repmat(experimentZero(m),nCon,1);
obj = objectiveZero([nCon,1]);
nTimes = length(times);
measurements = cell(nCon,1); % for plotting
for i = 1:nCon
    
    % Create sample experimental condition and generate some test data
    ici = [1, 2, 0]'; % same ics but generateTestData will randomize measuremens
    con_i = experimentInitialValue(m, ici, [], [], ['VariantCon' num2str(i)]); % modified ICs
    
    [outputsList, timesList, measurementsList] = generateTestData(m, con_i, times, outputs, sd);
    measurements{i} = reshape(measurementsList, nTimes, 3);
    
    obs = observationLinearWeightedSumOfSquares(outputsList, timesList, sd, ['VariantObs' num2str(i)]);
    obj_i = obs.Objective(measurementsList);
    
    newopts = [];
    newopts.UseParams = [i;1]; % different kf, same kr
    newopts.UseSeeds = [i;i;i]; % all different seeds
    newopts.ParamLowerBound = [1e-2 1e-2];
    newopts.ParamUpperBound = [1e2  1e2];
%     newopts.StartingParams = [4 4]; %%TODO%% implement this
    
    opts = BuildFitOpts(opts, m, con_i, obj_i, newopts);
    
    con(i) = con_i;
    obj(i) = obj_i;
    
end

% Plot test data
% figure
% for i = 1:nCon
%     subplot(1, nCon, i)
%     plot(times, measurements{i})
% end
% suptitle('Generated Test Data') % this function requires bioinformatics toolbox

% Fit multiple objective functions and experiments
opts.Verbose = 1;
opts.TolOptim = 1;
opts.MaxStepSize = 1;
opts.Normalized = false;
opts.UseAdjoint = true;
opts.RunParallelExpts = false;

[m, con, G, D] = FitObjective(m, con, obj, opts);
% end up with 3 models in m

% opts.RunParallelExpts = false;
% G = ObjectiveValue(fit, opts);
% opts.RunParallelExpts = true;
% G_ = ObjectiveValue(fit, opts);

% opts.RunParallelExpts = false;
% [D, G] = ObjectiveGradient(fit, opts);
% opts.RunParallelExpts = true;
% [D_, G_] = ObjectiveGradient(fit, opts);

% opts.RunParallelExpts = false;
% [H, D, G] = ObjectiveHessian(fit, opts);
% opts.RunParallelExpts = true;
% [H_, D_, G_] = ObjectiveHessian(fit, opts);

% opts.RunParallelExpts = false;
% [F, FAll] = ObjectiveInformation(fit, opts);
% opts.RunParallelExpts = true;
% [F_, FAll_] = ObjectiveInformation(fit, opts);

%% Display fit results
timesFine = linspace(0, tF, 100)';
simFits = [];
for i = 1:nCon
    simFit = SimulateSystem(m(i), con(i), tF);
    simFits = [simFits; simFit];
end

figure
hold on
% Plot test data
for i = 1:nCon
    plot(times, measurements{i}, '+')
    ax = gca;
    ax.ColorOrderIndex = 1;
end
% Plot fits
for i = 1:nCon
    plot(timesFine, simFits(i).y(timesFine)')
    ax = gca;
    ax.ColorOrderIndex = 1;
end
legend('A','B','C')
xlabel('Time')
ylabel('Amount')
title('Varying seeds, varying kf')

