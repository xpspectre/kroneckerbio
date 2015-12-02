function m = equilibrium_model_analytic
% Basic analytic model with seeds
m = InitializeModelAnalytic('Equilibrium');

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

m = AddReaction(m, 'r1', {'A', 'B'}, {'C'}, 'A*B*kf', 'C*kr');

m = FinalizeModel(m);