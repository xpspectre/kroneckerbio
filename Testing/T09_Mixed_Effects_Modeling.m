% Test/tutorial script for multiple participants/experiments with datasets/objective functions

clear; close all; clc
rng('default')

%% Construct equilibrium experiment A + B <-> C with seeds
m = InitializeModel('Equilibrium');

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

%% Multiple models and experiments (using old notation)
% Set up experiments and generate some test data to fit
% times = [0.1, 1];
tF = 1;
times = linspace(0, tF, 5)';
outputs = {'A','B','C'};
sd = sdLinear(0.05, 0.1);

nExpts = 4;
nTimes = length(times);
expts(nExpts,1) = experimentZero(m);
objs = objectiveZero([nExpts,nExpts]);
measurements = cell(nExpts,1); % for plotting
for i = 1:nExpts
%     ici = [1+rand, 1+rand, 0]'; % randomized ics
    ici = [1, 2, 0]'; % same ics but generateTestData will randomize measuremens
    expt = experimentInitialValue(m, [], ici, [], [], ['VariantExpt' num2str(i)]); % modified ICs
    expts(i) = expt;
    
    [outputsList, timesList, measurementsList] = generateTestData(m, expt, times, outputs, sd);
    measurements{i} = reshape(measurementsList, nTimes, 3);
    
    obs = observationLinearWeightedSumOfSquares(outputsList, timesList, sd, ['VariantObs' num2str(i)]);
    obj = obs.Objective(measurementsList);
    objs(i,i) = obj;
end

% Plot test data
figure
for i = 1:nExpts
    subplot(1, nExpts, i)
    plot(times, measurements{i})
end
suptitle('Generated Test Data') % this function requires bioinformatics toolbox

% Fit multiple objective functions and experiments
opts = [];
opts.Verbose = 2;
opts.TolOptim = 1;
opts.MaxStepSize = 1;
% opts.UseParams = true(nk, 1);
opts.UseParams = [1:nExpts;ones(1,nExpts)]; % more complicated parameter fitting
% opts.UseParams = [1:nExpts;1:nExpts]; % more complicated parameter fitting
opts.UseSeeds = true(ns, nExpts); % default is to fit all the seeds
opts.UseAdjoint = false; % adjoint fails with current setup
[mFit, exptFit] = FitObjective(m, expts, objs, opts);

%% Display fit results
timesFine = linspace(0, tF, 100)';
simFits = [];
for i = 1:nExpts
    simFit = SimulateSystem(mFit(i), exptFit(i), tF);
    simFits = [simFits; simFit];
end

figure
hold on
% Plot test data
for i = 1:nExpts
    plot(times, measurements{i}, '+')
    ax = gca;
    ax.ColorOrderIndex = 1;
end
% Plot fits
for i = 1:nExpts
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
