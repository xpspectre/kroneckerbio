function m = equilibrium_dimer_model
% Slightly more complicated mass action model with seeds
m = InitializeModelMassActionAmount('Equilibrium_Dimer');

m = AddCompartment(m, 'Solution', 3, 1);

m = AddSeed(m, 'sA',   2);
m = AddSeed(m, 'sB',   2);
m = AddSeed(m, 'sC',   0);
m = AddSeed(m, 'sA:A', 0);

m = AddState(m, 'A',   'Solution', 'sA');
m = AddState(m, 'B',   'Solution', 'sB');
m = AddState(m, 'C',   'Solution', 'sC');
m = AddState(m, 'A:A', 'Solution', 'sA:A');

m = AddOutput(m, 'A',   'A');
m = AddOutput(m, 'B',   'B');
m = AddOutput(m, 'C',   'C');
m = AddOutput(m, 'A:A', 'A:A');

m = AddParameter(m, 'kf', 5);
m = AddParameter(m, 'kr', 3);
m = AddParameter(m, 'kAf', 5);
m = AddParameter(m, 'kAr', 5);

m = AddReaction(m, 'r1', 'A', 'A', 'A:A', '', 'kAf', 'kAr');
m = AddReaction(m, 'r2', 'A:A', 'B', 'C', '', 'kf', 'kr');

m = FinalizeModel(m);