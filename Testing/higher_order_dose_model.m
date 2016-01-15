function [m, con, obj, opts] = higher_order_dose_model()
m = AnalyticModel();
m.AddCompartment('v', 3, 1);
m.AddState('x1', 'v', 's1^2');
m.AddState('x2', 'v', 's1*k3');
m.AddState('x3', 'v', 'k1^2');
m.AddState('x4', 'v', 's1*s2');
m.AddState('x5', 'v', 'k1*k2');
m.AddSeed('s1', 0);
m.AddSeed('s2', 0);
m.AddParameter('k1', 2);
m.AddParameter('k2', 3);
m.AddParameter('k3', 7);
m.AddStatesAsOutputs;
m.Finalize;

dos = DoseConstant(m, [5,8], 2:10);
con = experimentInitialValue(m, [], [], dos);

values = [
    5   0       6.1
    3   1       4.2
    3   1.1     4.2
    1   3.4     26
    1   10      200
    2   8.9     250
    4   7       200
    ];
obs = observationLinearWeightedSumOfSquares(values(:,1), values(:,2), sdLinear(0.5, 0.1));
obj = obs.Objective(values(:,3));

opts.Verbose = false;
opts.UseParams = 1:m.nk;
opts.UseSeeds = [];
opts.UseInputControls = [];
opts.UseDoseControls = 1:2;
