function m = equilibrium_dimer_model
% Slightly more complicated mass action model with seeds
m = InitializeModelMassActionAmount('Equilibrium_Dimer');

m = AddCompartment(m, 'Solution', 3, 1);

m = AddSeed(m, 'A_0',   2);
m = AddSeed(m, 'B_0',   2);
m = AddSeed(m, 'C_0',   0);
m = AddSeed(m, 'A:A_0', 0);

m = AddState(m, 'A',   'Solution', 'A_0');
m = AddState(m, 'B',   'Solution', 'B_0');
m = AddState(m, 'C',   'Solution', 'C_0');
m = AddState(m, 'A:A', 'Solution', 'A:A_0');

m = AddOutput(m, 'A',   'A');
m = AddOutput(m, 'B',   'B');
m = AddOutput(m, 'C',   'C');
m = AddOutput(m, 'A:A', 'A:A');

m = AddParameter(m, 'kf', 5);
m = AddParameter(m, 'kr', 3);
m = AddParameter(m, 'kAf', 5);
m = AddParameter(m, 'kAr', 5);

m = AddReaction(m, 'r1', {'A', 'A'}, {'A:A'}, 'kAf', 'kAr');
m = AddReaction(m, 'r2', {'A:A', 'B'}, {'C'}, 'kf', 'kr');

m = FinalizeModel(m);