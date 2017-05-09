function m = cascade_model(N)
% MAPK-like cascade model of variable size for benchmarking

m = InitializeModelAnalytic(sprintf('Cascade%i', N));

m = AddCompartment(m, 'v', 3, 1);

m = AddSeed(m, 'A0', 10);

m = AddState(m, 'A', 'v', 'A0'); % initial activator
for iN = 1:N
    m = AddState(m, sprintf('X%i', iN), 'v', 10);
    m = AddState(m, sprintf('Xp%i', iN), 'v', 0);
end

for iN = 1:N
    m = AddParameter(m, sprintf('kon%i', iN), 5);
    m = AddParameter(m, sprintf('koff%i', iN), 3);
end

m = AddReaction(m, {'X1 activation', 'X1 inactivation'}, {'X1'}, {'Xp1'}, 'kon1*X1*A', 'koff1*Xp1');
for iN = 2:N
    m = AddReaction(m ,{sprintf('X%i activation', iN), sprintf('X%i inactivation', iN)}, {sprintf('X%i', iN)}, {sprintf('Xp%i', iN)}, sprintf('kon%i*X%i*Xp%i', iN, iN, iN-1), sprintf('koff%i*Xp%i', iN, iN));
end

for iN = 1:N
    m = AddOutput(m, sprintf('X%i', iN), sprintf('X%i', iN));
    m = AddOutput(m, sprintf('Xp%i', iN), sprintf('Xp%i', iN));
end

m = FinalizeModel(m);
end

