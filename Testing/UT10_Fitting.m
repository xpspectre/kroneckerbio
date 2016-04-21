function tests = UT10_Fitting()
tests = functiontests(localfunctions);
if nargout < 1
    tests.run;
end
end

function testSimpleFitting(a)
[m, con, obj, opts] = simple_model();
opts.MaxIter = 2;

Gold = ObjectiveValue(m, con, obj, opts);

[m, con] = FitObjective(m, con, obj, opts);

Gnew = ObjectiveValue(m, con, obj, opts);

a.verifyLessThan(Gnew, Gold)
end

function testSimpleAnalyticFitting(a)
[m, con, obj, opts] = simple_analytic_model();
opts.MaxIter = 2;

Gold = ObjectiveValue(m, con, obj, opts);

[m, con] = FitObjective(m, con, obj, opts);

Gnew = ObjectiveValue(m, con, obj, opts);

a.verifyLessThan(Gnew, Gold)
end

function testMichaelisMentenFitting(a)
[m, con, obj, opts] = michaelis_menten_model();

Gold = ObjectiveValue(m, con, obj, opts);

timer = tic;
[m, con, G, D] = FitObjective(m, con, obj, opts);
time = toc(timer);

Gnew = ObjectiveValue(m, con, obj, opts);
Dnew = ObjectiveGradient(m, con, obj, opts);

a.verifyLessThan(Gnew, Gold)
a.verifyLessThan(time, 12) % sec, Hack to verify it finishes quickly with correct obj fun grad in FitObjective

a.verifyEqual(G, Gnew, 'RelTol', 1e-4)
a.verifyEqual(D, Dnew, 'RelTol', 1e-4)
end

function testMakeMultiConditionFit(a)
opts = [];
opts.Verbose = 0;
fit = FitObject('UT10_MultiConditionFit', opts);

m = equilibrium_model;
fit.addModel(m);

times = linspace(0, 1, 5)';
outputs = {'A','B','C'};
sd = sdLinear(0.05, 0.1);

nCon = 3;
for i = 1:nCon
    con = experimentInitialValue(m, [1,2,0]', [], [], ['Con' num2str(i)]);
    [outputsList, timesList, measurementsList] = generateTestData(m, con, times, outputs, sd);
    obs = observationLinearWeightedSumOfSquares(outputsList, timesList, sd, ['Obs' num2str(i)]);
    obj = obs.Objective(measurementsList);
    opts = [];
    opts.UseParams = [i;1];
    fit.addFitConditionData(obj, con, opts);
end

a.verifyEqual(length(fit.Models), nCon); % different UseParams -> added dummy models
a.verifyEqual(length(fit.Conditions), nCon);
a.verifyEqual(length(fit.Objectives), nCon);
end

function testMakeMultiModelFit(a)
opts = [];
opts.ModelsShareParams = true;
opts.Verbose = 0;

fit = FitObject('UT10_MultiModelFit', opts);

m1 = equilibrium_model;
fit.addModel(m1);
m2 = equilibrium_dimer_model;
fit.addModel(m2);

times = linspace(0, 1, 5)';
outputs = {'A','B','C'};
sd = sdLinear(0.05, 0.1);
    
con = experimentInitialValue(m1, [1,2,0], [], [], 'Con');
[outputsList, timesList, measurementsList] = generateTestData(m1, con, times, outputs, sd);
obs = observationLinearWeightedSumOfSquares(outputsList, timesList, sd, 'Obs1');
obj = obs.Objective(measurementsList);
opts = [];
opts.UseParams = [1,2]';
opts.UseSeeds = [1,2,0]';
fit.addFitConditionData(obj, con, opts);

con = experimentInitialValue(m2, [2,2,0,0], [], [], 'Con2');
[outputsList, timesList, measurementsList] = generateTestData(m2, con, times, outputs, sd);
obs = observationLinearWeightedSumOfSquares(outputsList, timesList, sd, 'Obs2');
obj = obs.Objective(measurementsList);
opts = [];
opts.UseParams = [3,2,4,5]';
opts.UseSeeds = [3,4,0,0]';
fit.addFitConditionData(obj, con, opts);

a.verifyEqual(length(fit.Models), 2);
a.verifyEqual(length(fit.Conditions), 2);
a.verifyEqual(length(fit.Objectives), 2);
end

% TODO: add FitObject-related tests
