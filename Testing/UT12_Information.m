function tests = UT12_Information()
tests = functiontests(localfunctions);
if nargout < 1
    tests.run;
end
end

function testObjectiveValueSimple(a)
[m, con, obj, opts] = simple_model();

F1 = ObjectiveInformation(m, con, obj, opts);
end

function testObjectiveInformationArguments(a)
% Test different ways of calling ObjectiveInformation
[m, con, obj, opts] = simple_model();

dxdTSol = {}; % blank for testing

[F1, FAll1] = ObjectiveInformation(m, con, obj, struct('Verbose', 0)); % only to suppress output; get FIM of default fit params
a.verifyEqual(size(F1), [10,10]);
a.verifyEqual(size(FAll1), [1,1]);

% Just see if the following 2 give the same result
[F2, FAll2] = ObjectiveInformation(m, con, obj, opts);
[F3, FAll3] = ObjectiveInformation(m, con, obj, opts, dxdTSol);
a.verifyEqual(F2, F3);
a.verifyEqual(FAll2, FAll3);
end

function testObjectiveValueSimpleAnalytic(a)
[m, con, obj, opts] = simple_analytic_model();

F1 = ObjectiveInformation(m, con, obj, opts);
end

function testObjectiveValueMichaelisMenten(a)
[m, con, obj, opts] = michaelis_menten_model();

F1 = ObjectiveInformation(m, con, obj, opts);
end