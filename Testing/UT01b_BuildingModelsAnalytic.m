function tests = UT01b_BuildingModelsAnalytic()
tests = functiontests(localfunctions);
if nargout < 1
    tests.run;
end
end

function testInitializingModel(a)
model_name = 'TestAnalytic';
m = AnalyticModel(model_name);
a.verifyEqual(m.Name, model_name);
end

function testInitializingModelWithEmptyName(a)
m = AnalyticModel();
a.verifyEqual(m.Name, '');
end

function testAddCompartment0(a)
m = AnalyticModel();
m.AddCompartment('test', 0, 1);
m.Finalize;
a.verifyEqual(m.Compartments(1).Dimension, 0);
end

function testNonUnityZeroCompartment(a)
m = AnalyticModel();
a.verifyError(@()m.AddCompartment('test0', 0, 2), 'KroneckerBio:Validate:CompartmentSize:NonUnityZeroCompartment');
end

function testAddCompartments(a)
m = AnalyticModel();
m.AddCompartment('test1', 1, 1);
m.AddCompartment('test2', 2, 1);
m.AddCompartment('test3', 3, 1);
m.Finalize;
a.verifyEqual(m.Compartments(1).Dimension, 1);
a.verifyEqual(m.Compartments(2).Dimension, 2);
a.verifyEqual(m.Compartments(3).Dimension, 3);
end

function testAddCompartment4(a)
m = AnalyticModel();
a.verifyError(@()m.AddCompartment('test4', 4, 1), 'KroneckerBio:Validator:CompartmentDimension');
end

function testStoichiometriyMatrixInRightOrder(a)
m = AnalyticModel();
m.AddCompartment('v1', 3, 1);
m.AddState('x1', 'v1', 2);
m.AddInput('u1', 'v1', 3);
m.AddParameter('k1', 4);
m.AddReaction('', 'x1', 'u1', 'x1*k1');
m.Finalize;

a.verifyEqual(m.m.f(0, 5, 3), -5*4)
end

function [m, x0, u0] = model_with_some_species()
m = AnalyticModel();
m.AddCompartment('v1', 3, 1);
m.AddState('x1', 'v1', 2);
m.AddState('x2', 'v1', 3);
m.AddState('x3', 'v1', 5);
m.AddState('x4', 'v1', 7);
m.AddParameter('k1', 4);
m.AddParameter('k2', 3);
m.Finalize;
x0 = m.m.x0(m.m.s);
u0 = vec([]);
end

function testAddReactionFwd(a)
[m, x0, u0] = model_with_some_species();
m.AddReaction('test', {'x1', 'x2'}, 'x3', 'k1*x1*x2'); % must be specified explicitly in analytic model
m.Finalize;
a.verifyEqual(m.m.f(0,x0,u0), [-24;-24;24;0])
end

function testAddReactionRev(a)
[m, x0, u0] = model_with_some_species();
m.AddReaction('test', {'x1', 'x2'}, 'x3', '', 'k2*x3');
m.Finalize;
a.verifyEqual(m.m.f(0,x0,u0), [15;15;-15;0])
end

function testRepeatedComponents(a)
% Consider adding repeated reaction and rule tests
m = AnalyticModel();
m.AddCompartment('v1', 3, 1);
m.AddParameter('k', 4);
m.AddSeed('s1', 1);
m.AddState('x1', 'v1', 1);
m.AddInput('u1', 'v1', 1);
m.AddOutput('o1', 'x1');

test = copy(m);
test.AddCompartment('v1', 3, 1);
a.verifyError(@()test.Finalize, 'KroneckerBio:FinalizeModel:RepeatCompartment');

test = copy(m);
test.AddParameter('k', 4);
a.verifyError(@()test.Finalize, 'KroneckerBio:FinalizeModel:RepeatParameter');

test = copy(m);
test.AddSeed('s1', 1);
a.verifyError(@()test.Finalize, 'KroneckerBio:FinalizeModel:RepeatSeed');

test = copy(m);
test.AddState('x1', 'v1', 1);
a.verifyError(@()test.Finalize, 'KroneckerBio:Parser:qualifyCompartments:AmbiguousSpeciesName'); % caught earlier than old RepeatSpecies exception

test = copy(m);
test.AddInput('u1', 'v1', 1);
a.verifyError(@()test.Finalize, 'KroneckerBio:FinalizeModel:RepeatSpecies');

test = copy(m);
test.AddOutput('o1', 'x1');
a.verifyError(@()test.Finalize, 'KroneckerBio:FinalizeModel:RepeatOutput');
end

function testDuplicateSpeciesNames(a)
m = AnalyticModel();
m.AddCompartment('v1', 3, 1);
m.AddCompartment('v2', 3, 1);
m.AddParameter('k', 4);
m.AddState('x1', 'v1');
m.AddState('x1', 'v2');
m.AddState('x2', 'v2');

a.verifyError(@()m.AddReaction('', {'', 'x1'}, {}, 'k*x1'), 'KroneckerBio:Validate:ReactionSpecies:InvalidBlankName')

a.verifyError(@()m.AddReaction('', {1, 'x1'}, {}, 'k*x1'), 'KroneckerBio:Validate:ReactionSpecies:InvalidName')

a.verifyError(@()m.AddReaction('', 'v1.x1.y1', {}, 'k*x1'), 'KroneckerBio:Validate:ReactionSpecies:InvalidNameDots')

a.verifyError(@()m.AddReaction('', 'v1.x1"y1', {}, 'k*x1'), 'KroneckerBio:Validate:ReactionSpecies:InvalidNameDoubleQuotes')

test = copy(m);
test.AddReaction('', 'v1.x1', {}, 'k*"v1.x1"');
test.Finalize;
a.verifyEqual(test.nr, 1) % no error on finalization

test = copy(m);
test.AddReaction('', 'x1', {}, 'k*x1', '', 'v1');
test.Finalize;
a.verifyEqual(test.nr, 1)

test = copy(m);
test.AddReaction('', {}, 'x1', 'k*x1', '', 'v1');
test.Finalize;
a.verifyEqual(test.nr, 1)

test = copy(m);
test.AddReaction('', 'x1', {}, 'k*"v1.x1"', '', 'v1');
test.Finalize;
a.verifyEqual(test.nr, 1)

test = copy(m);
test.AddReaction('', 'x1', {}, 'k*"v2.x2"', '', 'v1');
test.Finalize;
a.verifyEqual(test.nr, 1)

test = copy(m);
test.AddReaction('', 'v1.x2', {}, 'k*v1.x2');
a.verifyError(@()test.Finalize, 'KroneckerBio:Parser:qualifyCompartments:MissingQualifiedSpeciesName')

test = copy(m);
test.AddReaction('', 'x3', {}, 'k*x3');
a.verifyError(@()test.Finalize, 'KroneckerBio:Parser:qualifyCompartments:MissingUnqualifiedSpeciesName')

test = copy(m);
test.AddReaction('', 'x1', {}, 'k*x1');
a.verifyError(@()test.Finalize, 'KroneckerBio:Parser:qualifyCompartments:AmbiguousSpeciesName')

test = copy(m);
test.AddReaction('', 'x1', {}, 'k*"v1.x1"');
a.verifyError(@()test.Finalize, 'KroneckerBio:Parser:qualifyCompartments:AmbiguousSpeciesName')

test = copy(m);
test.AddReaction('', 'x1', {}, 'k*x2', '', 'v1');
a.verifyError(@()test.Finalize, 'KroneckerBio:Parser:qualifyCompartments:MissingSpeciesInReactionCompartment')

test = copy(m);
test.AddReaction('', 'x2', {}, 'k*"x2"', '', 'v1');
a.verifyError(@()test.Finalize, 'KroneckerBio:Parser:qualifyCompartments:MissingSpeciesInReactionCompartment')

test = copy(m);
test.AddReaction('', 'x1', {}, 'k*"x1"', '', 'v3');
a.verifyError(@()test.Finalize, 'KroneckerBio:AnalyticModel:qualifyNames:MissingReactionCompartment')
end

function testAddReactionsWithEmptyNames(a)
% Test various ways that reactions with empty reaction and product species lists
% are added. Does not test finalization.
m = AnalyticModel();
m.AddCompartment('v1', 3, 1);
m.AddParameter('k1', 4);
m.AddParameter('k2', 3);
m.AddState('x1', 'v1', 2);
m.AddState('x2', 'v1', 3);

%% Forms of empty species lists
emptySpeciesForms = {[], {}, '', cell(1,0)};
nForms = length(emptySpeciesForms);
nr = m.nr;

% Empty reactants
for i = 1:nForms
    m.AddReaction('r1', emptySpeciesForms{i}, 'x2', 'k1');
    a.verifyEqual(m.Reactions(nr+i).Reactants, cell(1,0))
end

% Empty products
for i = 1:nForms
    m.AddReaction('r1', 'x1', emptySpeciesForms{i}, 'k1');
    a.verifyEqual(m.Reactions(nr+nForms+i).Products, cell(1,0))
end

%% Forms of some empty species lists that should error
emptySpeciesForms = {{[]}, {''}};
nForms = length(emptySpeciesForms);

% Empty reactants
for i = 1:nForms
    a.verifyError(@()m.AddReaction('r1', emptySpeciesForms{i}, 'x2', 'k1'), 'KroneckerBio:Validate:ReactionSpecies:InvalidBlankName');
end

% Empty products
for i = 1:nForms
    a.verifyError(@()m.AddReaction('r1', 'x1', emptySpeciesForms{i}, 'k1'), 'KroneckerBio:Validate:ReactionSpecies:InvalidBlankName');
end
end

function testsymbolic2PseudoKroneckerMM(a)

m = michaelis_menten_model();

a.verifyEqual(m.nv, 1)
a.verifyEqual(m.nk, 2)
a.verifyEqual(m.ns, 1)
a.verifyEqual(m.nu, 1)
a.verifyEqual(m.nx, 2)
a.verifyEqual(m.ny, 3)

end

function testRemoveComponent(a)
m = AnalyticModel();
m.AddCompartment('v1', 3, 1);
m.AddParameter('k1', 4);
m.AddSeed('s1', 1);
m.AddState('x1', 'v1', 1);
m.AddInput('u1', 'v1', 1);
m.AddOutput('y1', 'x1');
m.AddReaction('r1', 'x1', {}, 'k1*x1');
m.Finalize;

test = copy(m);
test.RemoveCompartment('v1');
a.verifyEqual(test.nv, 0);
a.verifyEqual(length(test.Compartments), 0);

test = copy(m);
test.RemoveParameter('k1');
a.verifyEqual(test.nk, 0);
a.verifyEqual(length(test.Parameters), 0);

test = copy(m);
test.RemoveSeed('s1');
a.verifyEqual(test.ns, 0);
a.verifyEqual(length(test.Seeds), 0);

test = copy(m);
test.RemoveState('x1');
a.verifyEqual(test.nx, 0);
a.verifyEqual(length(test.States), 0);

test = copy(m);
test.RemoveInput('u1');
a.verifyEqual(test.nu, 0);
a.verifyEqual(length(test.Inputs), 0);

test = copy(m);
test.RemoveOutput('y1');
a.verifyEqual(test.ny, 0);
a.verifyEqual(length(test.Outputs), 0);

test = copy(m);
test.RemoveReaction('r1');
a.verifyEqual(test.nr, 0);
a.verifyEqual(length(test.Reactions), 0);

a.verifyError(@()m.RemoveCompartment('v2'), 'KroneckerBio:Model:RemoveCompartment:CompartmentNotFound');
a.verifyError(@()m.RemoveParameter('k2'), 'KroneckerBio:Model:RemoveParameter:ParameterNotFound');
a.verifyError(@()m.RemoveSeed('s2'), 'KroneckerBio:Model:RemoveSeed:SeedNotFound');
a.verifyError(@()m.RemoveState('x2'), 'KroneckerBio:Model:RemoveState:StateNotFound');
a.verifyError(@()m.RemoveInput('u2'), 'KroneckerBio:Model:RemoveInput:InputNotFound');
a.verifyError(@()m.RemoveOutput('y2'), 'KroneckerBio:Model:RemoveOutput:OutputNotFound');
a.verifyError(@()m.RemoveReaction('r2'), 'KroneckerBio:Model:RemoveReaction:ReactionNotFound');
end