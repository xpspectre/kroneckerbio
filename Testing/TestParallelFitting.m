% Test parallel fitting
% clear; close all; clc
% rng('default')

function TestParallelFitting
fitopts.MaxIter = 2;
nExperiments = 3;
nMeasurements = 15;
testfun = generateTestParallel('simple', fitopts, nExperiments, nMeasurements);
testfun();

end

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

    function test()
        % Tests whether parallel and serial FitObjective give the same
        % solution
        opts_parallel = opts;
        opts_parallel.RunParallelExpts = true;
        [mfit_parallel, confit_parallel, G_parallel, D_parallel] = ...
            FitObjective(m, repmat(con,nExperiments,1), obj, opts_parallel);

        [mfit, confit, G, D] = FitObjective(m, repmat(con,nExperiments,1), obj, opts);

%         testEqual = @(par, ser, message) a.verifyEqual(par, ser,...
%             'AbsTol', 1e-12, 'RelTol', 1e-9, message);

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