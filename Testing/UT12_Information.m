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
[m, con, obj, opts] = simple_model();

fit = FitObject.buildFitObject(m, con, obj, opts);

dxdTSol = {}; % blank for testing

[F1, All1] = ObjectiveInformation(fit);
a.verifyEqual(size(F1), [9,9]);
a.verifyEqual(size(All1), [1,1]);
% Assume the first 1 is correct. All the other invocations must return the exact same result

[F2, All2] = ObjectiveInformation(fit, opts);
a.verifyEqual(F2, F1);
a.verifyEqual(All2, All1);

[F3, All3] = ObjectiveInformation(fit, opts, dxdTSol);
a.verifyEqual(F3, F1);
a.verifyEqual(All3, All1);

% This one is excluded because opts.Verbose = 0 is required to suppress all
%   output and it species different params, seeds, etc. to fit
% [F4, All4] = ObjectiveInformation(m, con, obj);
% a.verifyEqual(F4, F1);

[F5, All5] = ObjectiveInformation(m, con, obj, opts);
a.verifyEqual(F5, F1);
a.verifyEqual(All5, All1);

[F6, All6] = ObjectiveInformation(m, con, obj, opts, dxdTSol);
a.verifyEqual(F6, F1);
a.verifyEqual(All6, All1);
end

function testObjectiveValueSimpleAnalytic(a)
[m, con, obj, opts] = simple_analytic_model();

F1 = ObjectiveInformation(m, con, obj, opts);
end

function testObjectiveValueMichaelisMenten(a)
[m, con, obj, opts] = michaelis_menten_model();

F1 = ObjectiveInformation(m, con, obj, opts);
end