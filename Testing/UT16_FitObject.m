function tests = UT16_FitObject()
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