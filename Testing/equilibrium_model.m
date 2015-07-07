function m = equilibrium_model
% Basic mass action model with seeds
m = InitializeModel('Equilibrium');

m = AddCompartment(m, 'Solution', 3, 1);

m = AddSeed(m, 'sA', 1);
m = AddSeed(m, 'sB', 2);
m = AddSeed(m, 'sC', 0);

m = AddState(m, 'A', 'Solution', 'sA');
m = AddState(m, 'B', 'Solution', 'sB');
m = AddState(m, 'C', 'Solution', 'sC');

m = AddOutput(m, 'A', 'A');
m = AddOutput(m, 'B', 'B');
m = AddOutput(m, 'C', 'C');

m = AddParameter(m, 'kf', 5);
m = AddParameter(m, 'kr', 3);

m = AddReaction(m, '', '', 'A', 'B', 'C', '', 'kf', 'kr');

m = FinalizeModel(m);