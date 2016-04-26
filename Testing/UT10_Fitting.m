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
opts = BuildFitOpts;
opts.Verbose = 0;

m = equilibrium_model;

times = linspace(0, 1, 5)';
outputs = {'A','B','C'};
sd = sdLinear(0.05, 0.1);

nCon = 3;
con = repmat(experimentZero(m),nCon,1);
obj = objectiveZero([nCon,1]);
for i = 1:nCon
    con_i = experimentInitialValue(m, [1,2,0]', [], [], ['Con' num2str(i)]);
    [outputsList, timesList, measurementsList] = generateTestData(m, con_i, times, outputs, sd);
    obs = observationLinearWeightedSumOfSquares(outputsList, timesList, sd, ['Obs' num2str(i)]);
    obj_i = obs.Objective(measurementsList);
    opts = BuildFitOpts(opts, m, con_i, obj_i, struct('UseParams', [i;1]));
end

a.verifyEqual(length(opts.fit.ModelNames), nCon);
a.verifyEqual(length(opts.fit.ConditionNames), nCon);
a.verifyEqual(length(opts.fit.ObjectiveNames), nCon);
a.verifyEqual(opts.fit.AddDummyModel, [0;1;1]); % different UseParams -> added dummy models
end

function testMakeMultiModelFit(a)
opts = BuildFitOpts;
opts.Verbose = 0;
opts.ModelsShareParams = true;

m1 = equilibrium_model;
m2 = equilibrium_dimer_model;

times = linspace(0, 1, 5)';
outputs = {'A','B','C'};
sd = sdLinear(0.05, 0.1);
    
con1 = experimentInitialValue(m1, [1,2,0], [], [], 'Con');
[outputsList, timesList, measurementsList] = generateTestData(m1, con1, times, outputs, sd);
obs = observationLinearWeightedSumOfSquares(outputsList, timesList, sd, 'Obs1');
obj1 = obs.Objective(measurementsList);
opts = BuildFitOpts(opts, m1, con1, obj1, struct('UseParams', [1;2], 'UseSeeds', [1;2;0]));

con2 = experimentInitialValue(m2, [2,2,0,0], [], [], 'Con2');
[outputsList, timesList, measurementsList] = generateTestData(m2, con2, times, outputs, sd);
obs = observationLinearWeightedSumOfSquares(outputsList, timesList, sd, 'Obs2');
obj2 = obs.Objective(measurementsList);
opts = BuildFitOpts(opts, m2, con2, obj2, struct('UseParams', [3;2;4;5], 'UseSeeds', [3;4;0;0]));

a.verifyEqual(length(opts.fit.ModelNames), 2);
a.verifyEqual(length(opts.fit.ConditionNames), 2);
a.verifyEqual(length(opts.fit.ObjectiveNames), 2);
end
