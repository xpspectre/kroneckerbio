% Quick and dirty script to run all the example scripts. Just make sure this
% runs w/o errors.

function tests = RunExampleTests()
tests = functiontests(localfunctions);
if nargout < 1
    tests.run;
end
end

function testT01a(a)
T01a_Building_Model_MassAction;
end

function testT01b(a)
T01b_Building_Model_Analytic;
end

function testT02(a)
T02_Modifying_Models;
end

function testT03a(a)
T03a_Loading_Models_Basic;
end

function testT03b(a)
T03b_Loading_Models_Advanced;
end

function testT04(a)
T04_Basic_ODE_Simulation;
end

function testT05a(a)
T05a_Fitting_Parameters_Basic;
end

function testT05b(a)
T05b_Mixed_Effects_Fitting;
end

function testT05c(a)
T05c_Fit_Multiple_Models;
end

function testT05d(a)
T05d_Nonlinear_Mixed_Effects_Fitting;
end

function testT05e(a)
% currently highly experimental/built-in SAEM broken
T05e_Builtin_NLME_Fitting;
end

function testT06(a)
% empty file
T06_Parameter_Uncertainty;
end

function testT07(a)
T07_Topology_Uncertainty;
end