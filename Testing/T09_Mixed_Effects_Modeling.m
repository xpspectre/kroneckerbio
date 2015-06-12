% Test/tutorial script for multiple participants/experiments with datasets/objective functions

clear; close all; clc
rng('default')

%% Initialize fit object
fit = FitObject('T09_Mixed_Effects_Fit');

%% Construct equilibrium experiment A + B <-> C with seeds
m = InitializeModel('Base');

m = AddCompartment(m, 'Solution', 3, 1);

m = AddSeed(m, 'sA', 1);
m = AddSeed(m, 'sB', 2);
m = AddSeed(m, 'sC', 0);

m = AddState(m, 'A', 'Solution', 'sA');
m = AddState(m, 'B', 'Solution', 'sB');
m = AddState(m, 'C', 'Solution', 'sC');

m = AddOutput(m, 'A', 'A');
m = AddOutput(m, 'B', 'B');
m = AddOutput(m, 'C', 'C');

m = AddParameter(m, 'kf', 5);
m = AddParameter(m, 'kr', 3);

m = AddReaction(m, '', '', 'A', 'B', 'C', '', 'kf', 'kr');

m = FinalizeModel(m);

ns = m.ns;
nk = m.nk;

opts = [];

fit.addModel(m, opts)

%% Set up experimental conditions and generate some test data for multi-"participant" fitting
tF = 1;
times = linspace(0, tF, 5)';
outputs = {'A','B','C'};
sd = sdLinear(0.05, 0.1);

nCon = 3;
nTimes = length(times);
measurements = cell(nCon,1); % for plotting
for i = 1:nCon
    
    % Create sample experimental condition and generate some test data
    ici = [1, 2, 0]'; % same ics but generateTestData will randomize measuremens
    con = experimentInitialValue(m, ici, [], [], ['VariantCon' num2str(i)]); % modified ICs
    
    [outputsList, timesList, measurementsList] = generateTestData(m, con, times, outputs, sd);
    measurements{i} = reshape(measurementsList, nTimes, 3);
    
    obs = observationLinearWeightedSumOfSquares(outputsList, timesList, sd, ['VariantObs' num2str(i)]);
    obj = obs.Objective(measurementsList);
    
    opts = [];
    opts.Verbose = 2;
    opts.UseParams = [i;1]; % different kf, same kr
    opts.UseSeeds = [i;i;i]; % all different seeds
    
    fit.addFitConditionData(obj, con, opts);
    
end

% Plot test data
figure
for i = 1:nCon
    subplot(1, nCon, i)
    plot(times, measurements{i})
end
suptitle('Generated Test Data') % this function requires bioinformatics toolbox

% Fit multiple objective functions and experiments
opts = [];
opts.Verbose = 2;
opts.TolOptim = 1;
opts.MaxStepSize = 1;
opts.UseAdjoint = false; % adjoint fails with current setup
fitOut = FitObjective(fit, opts);

%% Display fit results
timesFine = linspace(0, tF, 100)';
simFits = [];
for i = 1:nCon
    simFit = SimulateSystem(fitOut.Models(i), fitOut.Conditions(i), tF);
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
    plot(timesFine, simFits(i).y(timesFine,:)')
    ax = gca;
    ax.ColorOrderIndex = 1;
end
legend('A','B','C')
xlabel('Time')
ylabel('Amount')
title('Varying seeds, varying kf')

%% Construct model variants


%% Construct different subpopulation structure

%% Linearized parameter uncertainty
% F = ObjectiveInformation(mFit, con, obj);
