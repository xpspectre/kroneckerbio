function m = equilibrium_model_ode
% Basic analytic model using ODEs directly
m = InitializeModelAnalytic('EquilibriumODE');

m = AddCompartment(m, 'Solution', 3, 1);

m = AddSeed(m, 'A_0', 1);
m = AddSeed(m, 'B_0', 2);
m = AddSeed(m, 'C_0', 0);

m = AddState(m, 'A', 'Solution', 'A_0');
m = AddState(m, 'B', 'Solution', 'B_0');
m = AddState(m, 'C', 'Solution', 'C_0');

m = AddOutput(m, 'A', 'A');
m = AddOutput(m, 'B', 'B');
m = AddOutput(m, 'C', 'C');

m = AddParameter(m, 'kf', 5);
m = AddParameter(m, 'kr', 3);

% Equivalent reaction:
%   m = AddReaction(m, 'r1', {'A', 'B'}, {'C'}, 'kf', 'kr');
m = AddOde(m, 'A', '-kf*A*B + kr*C');
m = AddOde(m, 'B', '-kf*A*B + kr*C');
m = AddOde(m, 'C', 'kf*A*B - kr*C');

m = FinalizeModel(m);