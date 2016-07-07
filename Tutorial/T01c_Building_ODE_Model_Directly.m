%% T01c Building a Model Directly From ODEs
% Sometimes, it's useful or necessary to build up the ODE (right hand side)
%   directly instead of specifying the reactions.
% Everything is the same as building an analytic model, except you use
%   AddOde instead of AddReaction. The two cannot coexist in a single
%   model. Each state must have a corresponding ODE.
%
% This tutorial builds a basic model in 2 ways: using reactions and using
%   ODEs.
% More details for simulating models can be found in later tutorials.

%% Build base model components
m = InitializeModelAnalytic('Equilibrium');

m = AddCompartment(m, 'Solution', 3, 1);

m = AddSeed(m, 'A_0', 1);
m = AddSeed(m, 'B_0', 2);
m = AddSeed(m, 'C_0', 0);

m = AddState(m, 'A', 'Solution', 'A_0');
m = AddState(m, 'B', 'Solution', 'B_0');
m = AddState(m, 'C', 'Solution', 'C_0');

m = AddOutput(m, 'A', 'A');
m = AddOutput(m, 'B', 'B');
m = AddOutput(m, 'C', 'C');

m = AddParameter(m, 'kf', 5);
m = AddParameter(m, 'kr', 3);

%% Build analytic model with reactions
m1 = m;
m1 = AddReaction(m1, 'r1', {'A', 'B'}, {'C'}, 'A*B*kf', 'C*kr');
m1 = FinalizeModel(m1);

%% Build analytic model with ODEs
m2 = m;
m2 = AddOde(m2, 'A', '-kf*A*B + kr*C');
m2 = AddOde(m2, 'B', '-kf*A*B + kr*C');
m2 = AddOde(m2, 'C', 'kf*A*B - kr*C');
m2 = FinalizeModel(m2);

% Debug: look at a random derivative inside to verify
% m1.dfdx(1,[1;29;3],[])
% m2.dfdx(1,[1;29;3],[])

%% Simulate models
% They should give the same results
tF = 1; % final time

con1 = experimentInitialValue(m1, [], [], [], 'InitialValueExperiment');
sim1 = SimulateSystem(m1, con1, tF);
sim1_2 = SimulateSensitivity(m1, con1, tF);
sim1_3 = SimulateCurvature(m1, con1, tF);
sim1_4 = SimulateLna(m1, con1, observationAll(tF));

figure
times = linspace(0, tF, 100);
plot(times, sim1.x(times))
legend('A','B','C')
xlabel('Time')
ylabel('Amount')
title('Using Reactions Model')

con2 = experimentInitialValue(m2, [], [], [], 'InitialValueExperiment');
sim2 = SimulateSystem(m2, con2, tF);
sim2_2 = SimulateSensitivity(m2, con2, tF);
sim2_3 = SimulateCurvature(m2, con2, tF);
sim2_4 = SimulateLna(m2, con2, observationAll(tF));

figure
times = linspace(0, tF, 100);
plot(times, sim2.x(times))
legend('A','B','C')
xlabel('Time')
ylabel('Amount')
title('Using ODE Model')

%% Fit models
outputList   = [1    1    2    2    3    3  ]'; % index of the output this measurement refers to
timesList    = [0.1  1    0.1  1    0.1  1  ]';
measurements = [0.6  0.4  1.5  1.3  0.4  0.6]';
sd = sdLinear(0.05, 0.1);
obs = observationLinearWeightedSumOfSquares(outputList, timesList, sd, 'DefaultObservation');
obj = obs.Objective(measurements);

%% Fit reactions
mFit1 = FitObjective(m1, con1, obj);
tF = 1;
times = linspace(0, tF, 100);
simOriginal1 = SimulateSystem(m1, con1, tF);
simFit1 = SimulateSystem(mFit1, con1, tF);

figure
hold on
plot(times, simOriginal1.x(times))
ax = gca;
ax.ColorOrderIndex = 1;
plot(timesList(1:2), reshape(measurements,2,3), '+')
ax = gca;
ax.ColorOrderIndex = 1;
plot(times, simFit1.x(times), ':')
hold off
legend('A','B','C','A data','B data','C data','A fit','B fit','C Fit')
xlabel('Time')
ylabel('Amount')
title('Using Reactions Model')

%% Fit ODEs
mFit2 = FitObjective(m2, con2, obj);
tF = 1;
times = linspace(0, tF, 100);
simOriginal2 = SimulateSystem(m2, con2, tF);
simFit2 = SimulateSystem(mFit2, con2, tF);

figure
hold on
plot(times, simOriginal2.x(times))
ax = gca;
ax.ColorOrderIndex = 1;
plot(timesList(1:2), reshape(measurements,2,3), '+')
ax = gca;
ax.ColorOrderIndex = 1;
plot(times, simFit2.x(times), ':')
hold off
legend('A','B','C','A data','B data','C data','A fit','B fit','C Fit')
xlabel('Time')
ylabel('Amount')
title('Using ODE Model')
