function tests = UT19_FitObject()
tests = functiontests(localfunctions);
if nargout < 1
    tests.run;
end
end

%% File fixture functions
function setupOnce(a)
m = equilibrium_model;
% m = equilibrium_model_analytic;
a.TestData.model = m;
a.TestData.experiment = experimentInitialValue(m, [], [], [], 'InitialValueExperiment');

outputList   = [1    1    2    2    3    3  ]'; % index of the output this measurement refers to
timesList    = [0.1  1    0.1  1    0.1  1  ]';
measurements = [0.6  0.4  1.5  1.3  0.4  0.6]';
sd = sdLinear(0.05, 0.1); % Simple error model with a constant and proportional term
obs = observationLinearWeightedSumOfSquares(outputList, timesList, sd, 'DefaultObservation');
a.TestData.objfun = obs.Objective(measurements);

a.TestData.experiment2 = experimentInitialValue(m, [], [], [], 'InitialValueExperiment2');

outputList   = [1    1    2    2    3    3  ]';
timesList    = [0.1  1    0.1  1    0.1  1  ]';
measurements = [0.7  0.5  1.6  1.4  0.5  0.7]';
sd = sdLinear(0.05, 0.1);
obs = observationLinearWeightedSumOfSquares(outputList, timesList, sd, 'DefaultObservation2');
a.TestData.objfun2 = obs.Objective(measurements);

end

function teardownOnce(a)
% pass
end

%% Tests
function testInitFitObject(a)
name = 'TestFit';
fit = FitObject(name);
a.assertMatches(fit.Name, name);
end

function testBuildFitObject(a)
fit = FitObject('TestFit');

fit.addModel(a.TestData.model);
a.assertEqual(size(fit.Models), [1,1])
a.assertMatches(fit.Models.Name, a.TestData.model.Name);

fit.addFitConditionData(a.TestData.objfun, a.TestData.experiment);
a.assertEqual(size(fit.Conditions), [1,1])
a.assertEqual(size(fit.Objectives), [1,1])
a.assertMatches(fit.Conditions.Name, a.TestData.experiment.Name);
a.assertMatches(fit.Objectives.Name, a.TestData.objfun.Name);

% Test parameter collecting
[T, details] = fit.collectParams;
a.assertEqual(numel(T), 2);
a.assertEqual(numel(T), height(details));
a.assertTrue(isequal(details{:,1}, {'kf';'kr'}));
end

function testParamSpec(a)
% Test easier-to-use API that lets you specify parameters to fit, their bounds,
% and how they should be shared between individuals
fit = FitObject('TestFit');
fit.addModel(a.TestData.model);

% Test that 1st fit condition data works
Name   = {'kf', 'kr'}';
Type   = {'k',  'k'}';
LB = 10.^[-2,   -2]';
UB = 10.^[ 6,    6]';
Use    = [ 1,    0]';
opts = [];
opts.ParamSpec = table(Type, LB, UB, Use, 'RowNames', Name); % variable names are necessary - they get turned into column names
fit.addFitConditionData(a.TestData.objfun, a.TestData.experiment, opts);
a.assertEqual(fit.componentMap, [1,1,1]);

% Test that additional fit condition data that needs to create a dummy model works
Use    = [ 2,    0]';
opts = [];
opts.ParamSpec = table(Type, LB, UB, Use, 'RowNames', Name);
fit.addFitConditionData(a.TestData.objfun2, a.TestData.experiment2, opts);
a.assertEqual(fit.componentMap, [1,1,1;2,2,2]);

end

function testMassActionModelUpdateExtra(a)
% Make sure extra field works and updates persist across m.Update
m = a.TestData.model;

[extra, extra2] = makeExtraData;

% Extra data survives m.Update
m = m.UpdateExtra(extra);
k = m.k;
knew = k + rand(size(k));
m = m.Update(knew);
verifyExtraData(a, m.Extra, extra);

% Extra data can be overridden
m = m.UpdateExtra(extra2);
verifyExtraData(a, m.Extra, extra, extra2);
end

function testAnalyticModelUpdateExtra(a)
m = equilibrium_model_analytic;

[extra, extra2] = makeExtraData;

% Extra data survives m.Update
m = m.UpdateExtra(extra);
k = m.k;
knew = k + rand(size(k));
m = m.Update(knew);
verifyExtraData(a, m.Extra, extra);

% Extra data can be overridden
m = m.UpdateExtra(extra2);
verifyExtraData(a, m.Extra, extra, extra2);
end

function testInitialValueExperimentUpdateExtra(a)
con = a.TestData.experiment;

[extra, extra2] = makeExtraData;

% Extra data survives con.Update
con = con.UpdateExtra(extra);
s = con.s;
snew = s + rand(size(s));
con = con.Update(snew, [], []);
verifyExtraData(a, con.Extra, extra);

% Extra data can be overridden
con = con.UpdateExtra(extra2);
verifyExtraData(a, con.Extra, extra, extra2);
end

function testSteadyStateExperimentUpdateExtra(a)
m = a.TestData.model;
con = experimentSteadyState(m, [], [], [], [], [], 'SteadyStateExperiment');

[extra, extra2] = makeExtraData;

% Extra data survives con.Update
con = con.UpdateExtra(extra);
s = con.s;
snew = s + rand(size(s));
con = con.Update(snew, [], []);
verifyExtraData(a, con.Extra, extra);

% Extra data can be overridden
con = con.UpdateExtra(extra2);
verifyExtraData(a, con.Extra, extra, extra2);
end

function [extra, extra2] = makeExtraData
extra = [];
extra.a = 1;
extra.b = [1,2,3];
extra.c = {1,2,3};
extra.d = 'extrastring';

extra2 = [];
extra2.b = [9,9,9];
end

function verifyExtraData(a, teststruct, extra, extra2)
if nargin == 4
    extra = pastestruct(extra, extra2);
end
teststruct = pastestruct(extra, teststruct); % just look for fields matching extra
a.verifyEqual(teststruct, extra);
end