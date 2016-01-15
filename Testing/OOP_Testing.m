% Test OOP model building
clear; close all; clc
rng('default');

%% Massaction model
m1 = MassActionAmountModel('test_massaction');

m1.AddCompartment('v1', 3, 1);
m1.AddCompartment('v2', 3, 1);

m1.AddSeed('s1', 5);

m1.AddState('x1', 'v1');
m1.AddState('x2', 'v1', 10);
m1.AddState('x3', 'v1', 0);
m1.AddState('x4', 'v1', 's1');

m1.AddInput('u1', 'v1', 1);
m1.AddInput('u2', 'v1', 2);
m1.AddInput('u3', 'v1', 3);
m1.AddInput('u4', 'v1', 4);

m1.AddOutput('y1', 'x1');
m1.AddOutput('y2', {'x1', 2});

m1.AddParameter('k1', 1);
m1.AddParameter('k2', 2);
m1.AddParameter('k3', 3);

m1.AddReaction('r1', {}, {}, 'k1');
m1.AddReaction('r2', {}, 'x3', 'k1');
m1.AddReaction('r3', {'x1'}, {'x3', 'x4'}, 'k1');

m1.Finalize;

% Copy massaction model
m1b = copy(m1);

% Modify models and make sure they remain independent
kNew = [1.1 2.2 3.3]';
m1.Update(kNew);

[m1.Parameters.Value]'
[m1b.Parameters.Value]'

%% Analytic model
m2 = AnalyticModel('test_analytic');

m2.AddCompartment('v1', 3, 1);
m2.AddCompartment('v2', 3, 1);

m2.AddSeed('s1', 5);

m2.AddState('x1', 'v1');
m2.AddState('x2', 'v1', 10);
m2.AddState('x3', 'v1', 0);
m2.AddState('x4', 'v1', 's1');

m2.AddInput('u1', 'v1', 1);
m2.AddInput('u2', 'v1', 2);
m2.AddInput('u3', 'v1', 3);
m2.AddInput('u4', 'v1', 4);

m2.AddOutput('y1', 'x1');
m2.AddOutput('y2', '2*x1');
m2.AddOutput('y3', 'k1*x1');

m2.AddParameter('k1', 1);
m2.AddParameter('k2', 2);
m2.AddParameter('k3', 3);

m2.AddReaction('r1', {}, {}, 'k1');
m2.AddReaction('r2', {}, 'x3', 'k1');
m2.AddReaction('r3', {'x1'}, {'x3', 'x4'}, 'k1*x1');

m2.AddRule('k4', '2*k2', 'repeated assignment'); % different from SBML rules in that they can't override/define exsiting params, can only make new ones; isn't really compatible with defining rules in terms of rules

m2.Finalize;

% Copy analytic model
m2b = copy(m2);

% Modify models and make sure they remain independent
kNew = [1.2 2.3 3.4]';
m2.Update(kNew);

[m2.Parameters.Value]'
[m2b.Parameters.Value]'

%% General actions
m = [m1,m2]; % compose different model types into a matrix