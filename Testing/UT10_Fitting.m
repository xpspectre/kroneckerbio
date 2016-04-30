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

function testSimpleParallelFitting(a)

fitopts.MaxIter = 2;
nExperiments = 3;
nTotalTimePoints = 15;
testfun = generateTestParallel('simple', fitopts, nExperiments, nTotalTimePoints);
testfun(a)
end

function testMichaelisMentenParallelFitting(a)

fitopts.MaxIter = 2;
nExperiments = 3;
nTotalTimePoints = 15;
testfun = generateTestParallel('michaelis_menten', fitopts, nExperiments, nTotalTimePoints);
testfun(a)
end

%% Test generation function for parallel tests
function testfun = generateTestParallel(model, fitopts, nExperiments, nMeasurements)

% Check for parallel toolbox. If missing, just skip this test
noParallelToolbox = isempty(ver('distcomp'));
if noParallelToolbox
    testfun = @(a)true;
    return
end

switch model
    case 'simple'
        [m, con, obj, opts] = simple_model();
    case 'simple_analytic'
        [m, con, obj, opts] = simple_analytic_model();
    case 'michaelis_menten'
        [m, con, obj, opts] = michaelis_menten_model();
end

% Get randomized experiments and their fitting data
[~, ~, obj, opts] = randomExperimentFittingData(m, con, opts, ...
    nExperiments, obj.tF, nMeasurements);

opts = mergestruct(opts, fitopts);

testfun = @test;

    function test(a)
        % Tests whether parallel and serial FitObjective give the same
        % solution
        opts_parallel = opts;
        opts_parallel.RunParallelExpts = true;
        [mfit_parallel, confit_parallel, G_parallel, D_parallel] = ...
            FitObjective(m, repmat(con,nExperiments,1), obj, opts_parallel);

        [mfit, confit, G, D] = FitObjective(m, repmat(con,nExperiments,1), obj, opts);

        testEqual = @(par, ser, message) a.verifyEqual(par, ser,...
            'AbsTol', 1e-12, 'RelTol', 1e-9, message);

        testEqual(G_parallel, G, 'Objective function values were not equal')
        testEqual(D_parallel, D, 'Gradients were not equal')
        testEqual(mfit_parallel.k, mfit.k, 'k values were not equal')
        for i_con = 1:nExperiments
            testEqual(confit_parallel(i_con).s, confit(i_con).s, ['s values were not equal in experiment ' num2str(i_con)])
            testEqual(confit_parallel(i_con).q, confit(i_con).q, ['q values were not equal in experiment ' num2str(i_con)])
            testEqual(confit_parallel(i_con).h, confit(i_con).h, ['h values were not equal in experiment ' num2str(i_con)])
        end
    end

end
