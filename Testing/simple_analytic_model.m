function [m, con, obj, opts] = simple_analytic_model()

m = LoadModelSbmlAnalytic('simple_analytic.xml');

% Convert kinetic parameter to seed
ind = lookupmember('growth_factor_injection', {m.Parameters(1:m.nk).Name});
m = AddSeed(m, 'growth_factor_injection', 8);
m = RemoveParameter(m, 'growth_factor_injection');
m = addStatesAsOutputs(m);
m = FinalizeModel(m);

dos = DoseConstant(m, 8, [10,20]);
inp = inputConstant(m, 2.5);
con = experimentInitialValue(m, [], inp, dos);

values = [
    4,  5,  0.25
    4,  40, 0.35
    5,  7,  2.1
    5, 25,  6.8
    7,  1,  9.1
    7, 30,  25
    ];

sd = sdLinear(0.5, 0.1);
obs = observationLinearWeightedSumOfSquares(values(:,1), values(:,2), sd);
obj = obs.Objective(values(:,3));

opts.Verbose = false;
opts.UseParams = 1:2:12;
opts.UseSeeds = 1;
opts.UseInputControls = [];
opts.UseDoseControls = [];
