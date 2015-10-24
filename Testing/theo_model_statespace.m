function m = theo_model_statespace
% Theophylline oral dose model expressed directly using ODE RHSs

m = InitializeModelAnalytic('TheoModelStateSpace');

m = AddCompartment(m, 'v', 3, 1);

m = AddSeed(m, 'Xg0', 320); % mg, bolus dose at t = 0 represented as a an initial condition/seed of GI compartment
m = AddSeed(m, 'Xp0', 0); % mg

m = AddState(m, 'Xg', 'v', 'Xg0'); % mg
m = AddState(m, 'Xp', 'v', 'Xp0'); % mg

m = AddParameter(m, 'ka',  1.5); % 1/hr
m = AddParameter(m, 'kel', 0.15); % 1/hr
m = AddParameter(m, 'V',   25); % L

m = AddStateSpaceEqAnalytic(m, 'Xg', '-ka*Xg');
m = AddStateSpaceEqAnalytic(m, 'Xp', 'ka*Xg - kel*Xp');

m = AddOutput(m, 'Cp', 'Xp/V'); % mg/L

m = FinalizeModel(m);

end