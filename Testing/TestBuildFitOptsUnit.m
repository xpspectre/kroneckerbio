clear; close all; clc
rng('default');

%% Test model addition
opts = [];

% Add 1st con
m = struct('Name', 'TestModel1', 'nk', 2);
con = struct('Name', 'TestCondition1', 'ns', 3, 'nq', 1, 'nh', 1);
obj = struct('Name', 'TestObjective1');
newopts = [];

opts = BuildFitOpts(opts, m, con, obj, newopts);

% Add 2nd con with different model
m2 = struct('Name', 'TestModel2', 'nk', 3, 'ns', 4);
con2 = struct('Name', 'TestCondition2', 'ns', 3, 'nq', 1, 'nh', 1);
obj2 = struct('Name', 'TestObjective2');
newopts2 = [];

opts = BuildFitOpts(opts, m2, con2, obj2, newopts2);

% Add 3rd con with same model as con1 and same UseParams
con3 = struct('Name', 'TestCondition3', 'ns', 3, 'nq', 1, 'nh', 1);
obj3 = struct('Name', 'TestObjective3');
newopts3 = [];

opts = BuildFitOpts(opts, m, con3, obj3, newopts3);

% Add 4th con with same model as con1 and different UseParams
con4 = struct('Name', 'TestCondition4', 'ns', 3, 'nq', 1, 'nh', 1);
obj4 = struct('Name', 'TestObjective4');
newopts4 = [];
newopts4.UseParams = [1;2]; % 2nd param is independent

opts = BuildFitOpts(opts, m, con4, obj4, newopts4);

% Add 5th con with same model as con1 and different UseParams and 2 obj funs
con5 = struct('Name', 'TestCondition5', 'ns', 3, 'nq', 1, 'nh', 1);
obj5_1 = struct('Name', 'TestObjective5_1');
obj5_2 = struct('Name', 'TestObjective5_2');
newopts5 = [];
newopts5.UseParams = [1;3]; % 2nd param is independent

opts = BuildFitOpts(opts, m, con5, [obj5_1; obj5_2], newopts5);

1;