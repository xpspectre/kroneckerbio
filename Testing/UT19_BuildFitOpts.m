function tests = UT19_BuildFitOpts()
tests = functiontests(localfunctions);
if nargout < 1
    tests.run;
end
end

%% File fixture functions
function setupOnce(a)
m = equilibrium_model;
m2 = equilibrium_model_analytic;
a.TestData.m = m;
a.TestData.m2 = m2;

a.TestData.con1 = experimentInitialValue(m, [], [], [], 'Con1');

outputList   = [1    1    2    2    3    3  ]'; % index of the output this measurement refers to
timesList    = [0.1  1    0.1  1    0.1  1  ]';
measurements = [0.6  0.4  1.5  1.3  0.4  0.6]';
sd = sdLinear(0.05, 0.1); % Simple error model with a constant and proportional term
obs = observationLinearWeightedSumOfSquares(outputList, timesList, sd, 'Obs1');
a.TestData.obj1 = obs.Objective(measurements);

a.TestData.con2 = experimentInitialValue(m, [], [], [], 'Con2');

outputList   = [1    1    2    2    3    3  ]';
timesList    = [0.1  1    0.1  1    0.1  1  ]';
measurements = [0.7  0.5  1.6  1.4  0.5  0.7]';
sd = sdLinear(0.05, 0.1);
obs = observationLinearWeightedSumOfSquares(outputList, timesList, sd, 'Obs2');
a.TestData.obj2 = obs.Objective(measurements);

end

function teardownOnce(a)
% pass
end

%% Tests
function testInitBuildFitOpts(a)
opts = BuildFitOpts;
a.assertTrue(isfield(opts, 'fit'));
end

function testBuildFitObject(a)
opts = BuildFitOpts;
opts = BuildFitOpts(opts, a.TestData.m, a.TestData.con1, a.TestData.obj1, []);

a.assertEqual(size(opts.fit.ModelNames), [1,1])
a.assertEqual(size(opts.fit.UseParams), [1,1])
a.assertEqual(opts.fit.UseParams{1}, [1;1]);

a.assertEqual(size(opts.fit.ConditionNames), [1,1])
a.assertEqual(size(opts.fit.ObjectiveNames), [1,1])
a.assertEqual(size(opts.fit.UseSeeds), [1,1])
a.assertEqual(opts.fit.UseSeeds{1}, [0;0;0])

a.assertEqual(opts.fit.ComponentMap, {1,1,1})

% Test parameter collecting - can't do this now because it's private
% [T, details] = collectAllActiveParameters(a.TestData.m, a.TestData.con1, opts.fit.ComponentMap, opts.fit.Tlocal2T);
% a.assertEqual(numel(T), 2);
% a.assertEqual(numel(T), height(details));
% a.assertTrue(isequal(details{:,1}, {'kf';'kr'}));
end

function testParamSpec(a)
% Test easier-to-use API that lets you specify parameters to fit, their bounds,
% and how they should be shared between individuals
opts = BuildFitOpts;

% Test that 1st fit condition data works
Name   = {'kf', 'kr'}';
Type   = {'k',  'k'}';
LB = 10.^[-2,   -2]';
UB = 10.^[ 6,    6]';
Use    = [ 1,    0]';
opts = BuildFitOpts(opts, a.TestData.m, a.TestData.con1, a.TestData.obj1, struct('ParamSpec', table(Type, LB, UB, Use, 'RowNames', Name)));
a.assertEqual(opts.fit.ComponentMap, {1,1,1})

% Test that additional fit condition data that needs to create a dummy model works
Use    = [ 2,    0]';
opts = BuildFitOpts(opts, a.TestData.m, a.TestData.con2, a.TestData.obj2, struct('ParamSpec', table(Type, LB, UB, Use, 'RowNames', Name)));
a.assertEqual(opts.fit.ComponentMap, {1,1,1;1,2,2});
a.assertEqual(opts.fit.AddDummyModel, [0;1]);
end

function testMassActionModelUpdateField(a)
% Make sure UpdateField works
m = a.TestData.m;

% Extra data survives m.Update
m = m.UpdateField(struct('Name', 'ModelNewName'));
k = m.k;
knew = k + rand(size(k));
m = m.Update(knew);
a.assertMatches(m.Name, 'ModelNewName');

% Extra data can be overridden
m = m.UpdateField(struct('Name', 'ModelNewName2'));
a.assertMatches(m.Name, 'ModelNewName2');
end

function testAnalyticModelUpdateField(a)
m = a.TestData.m2;

% Extra data survives m.Update
m = m.UpdateField(struct('Name', 'ModelNewName'));
k = m.k;
knew = k + rand(size(k));
m = m.Update(knew);
a.assertMatches(m.Name, 'ModelNewName');

% Extra data can be overridden
m = m.UpdateField(struct('Name', 'ModelNewName2'));
a.assertMatches(m.Name, 'ModelNewName2');
end

function testSingleModelTypeFit(a)
% Partially redundant with testMakeMultiConditionFit in UT10
opts = BuildFitOpts();

% Add 1st con
m = a.TestData.m;
con1 = experimentInitialValue(m, [], [], [], 'InitialValueExperiment1');
outputList   = [1    1    2    2    3    3  ]'; % index of the output this measurement refers to
timesList    = [0.1  1    0.1  1    0.1  1  ]';
measurements = [0.6  0.4  1.5  1.3  0.4  0.6]';
sd = sdLinear(0.05, 0.1);
obs = observationLinearWeightedSumOfSquares(outputList, timesList, sd, 'DefaultObservation1');
obj1 = obs.Objective(measurements);

newopts = [];

opts = BuildFitOpts(opts, m, con1, obj1, newopts);

% Add 2nd con with same model as con1 and same UseParams, different seeds
con2 = experimentInitialValue(m, [1.5;1.5;0], [], [], 'InitialValueExperiment2');
outputList   = [1    1    2    2    3    3  ]'; % index of the output this measurement refers to
timesList    = [0.1  1    0.1  1    0.1  1  ]';
measurements = [0.6  0.4  1.5  1.3  0.4  0.6]';
sd = sdLinear(0.05, 0.1);
obs = observationLinearWeightedSumOfSquares(outputList, timesList, sd, 'DefaultObservation2');
obj2 = obs.Objective(measurements);
newopts2 = [];

opts = BuildFitOpts(opts, m, con2, obj2, newopts2);

% Add 3rd con with same model as con1 and different UseParams -> add dummy model
con3 = experimentInitialValue(m, [1.6;1.6;0], [], [], 'InitialValueExperiment3');
outputList   = [1    1    2    2    3    3  ]'; % index of the output this measurement refers to
timesList    = [0.1  1    0.1  1    0.1  1  ]';
measurements = [0.6  0.4  1.5  1.3  0.4  0.6]';
sd = sdLinear(0.05, 0.1);
obs = observationLinearWeightedSumOfSquares(outputList, timesList, sd, 'DefaultObservation3');
obj3 = obs.Objective(measurements);
newopts3 = [];
newopts3.UseParams = [1;2]; % 2nd param is independent

opts = BuildFitOpts(opts, m, con3, obj3, newopts3);

% Add 4th con with same model as con1 and different UseParams and 2 obj funs ->
%   add dummy model, 2 rows in opts.fit.ComponentMap
con4 = experimentInitialValue(m, [1.7;1.7;0], [], [], 'InitialValueExperiment4');
outputList   = [1    1    2    2    3    3  ]'; % index of the output this measurement refers to
timesList    = [0.1  1    0.1  1    0.1  1  ]';
measurements = [0.6  0.4  1.5  1.3  0.4  0.6]';
sd = sdLinear(0.05, 0.1);
obs = observationLinearWeightedSumOfSquares(outputList, timesList, sd, 'DefaultObservation4_1');
obj4_1 = obs.Objective(measurements);
measurements = [0.6  0.4  1.5  1.3  0.4  0.7]';
obs = observationLinearWeightedSumOfSquares(outputList, timesList, sd, 'DefaultObservation4_1');
obj4_2 = obs.Objective(measurements);
newopts4 = [];
newopts4.UseParams = [1;3]; % 2nd param is independent

opts = BuildFitOpts(opts, m, con4, [obj4_1; obj4_2], newopts4);

% Test that opts' fields are right
a.assertEqual(opts.fit.ComponentMap, {1,1,1;1,2,2;1,3,3;1,4,[4;5]});
a.assertEqual(opts.fit.AddDummyModel, [0;0;1;1]);

% Test getting and setting params, bounds, etc
opts.fit.T2Tlocal([4;5;6;7]);

% Test running the fit
% con = [con1; con2; con3; con4];
% obj = [obj1; obj2; obj3; obj4_1; obj4_2];
% fit = FitObjective(m, con, obj, opts); % this gives an error: computeObjSendAdj (line 119) Subscript indices must either be real positive integers or logicals.
end

function testSplitComponentMap(a)
cmap = [1,1,1;2,2,2;3,3,3]; % more workers than conditions
cmaps = splitComponentMap(cmap,6);
a.assertEqual(cmaps, {[1,1,1],[2,2,2],[3,3,3],zeros(0,3),zeros(0,3),zeros(0,3)}');

cmap = [1,1,1;2,2,2;3,3,3]; % more conditions than workers
cmaps = splitComponentMap(cmap,2);
a.assertEqual(cmaps, {[1,1,1;2,2,2],[3,3,3]}');

cmap = [1,1,1;2,2,2;3,3,3;3,3,4;3,3,5]; % multiple objectives per condition with more workers than conditions
cmaps = splitComponentMap(cmap,6);
a.assertEqual(cmaps, {[1,1,1],[2,2,2],[3,3,3;3,3,4;3,3,5],zeros(0,3),zeros(0,3),zeros(0,3)}');

cmap = [1,1,1;2,2,2;3,3,3;3,3,4;3,3,5]; % multiple objectives per condition with more conditions than workers
cmaps = splitComponentMap(cmap,2);
a.assertEqual(cmaps, {[1,1,1;2,2,2],[3,3,3;3,3,4;3,3,5]}');
end