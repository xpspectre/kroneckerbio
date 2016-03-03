function tests = UT17_Nlme()
tests = functiontests(localfunctions);
if nargout < 1
    tests.run;
end
end

%% File fixture functions
function setupOnce(a)

% Build sample NLME model
m = InitializeModelAnalytic('LinearModel');
m = AddCompartment(m, 'v', 3, 1);
m = AddSeed(m, 'x0', 1);
m = AddState(m, 'x', 'v', 'x0');
m = AddParameter(m, 'k', 1.5);
m = AddReaction(m, 'Production', {}, {'x'}, 'k');
m = AddOutput(m, 'y', 'x');

Omega = 0.2; % Omega specified as omegas^2 (for now)
Sigma = sqrt(0.05); % Sigma specified as parameters in R

m = AddEta(m, 'k');
m = AddOmega(m, {'k'}, Omega);

m = AddErrorModel(m, {'y'}, {'sigma__y^2'}, {'sigma__y'}, Sigma);

m = FinalizeModel(m);

a.TestData.model = m;

%TODO: Make experiment and obj fun setup components

end

function teardownOnce(a)
% pass
end

%% Tests
function testBuildNlmeModel(a)
% Note that the first actual test in this file actually tests if the NLME model
% in the setupOnce function builds. So if it doesn't, this will fail as well.
a.assertTrue(true);
end

function testNlmeModelHasComponents(a)
m = a.TestData.model;

a.assertTrue(isfield(m, 'OmegaNames'));
a.assertTrue(isfield(m, 'fOmega'));

k = m.k;
k = rand(size(k));
m = m.Update(k);

a.assertTrue(isfield(m, 'fOmega')); % Makes sure that fOmega is in the model, that it isn't removed by m.Update

end