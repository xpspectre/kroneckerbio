function m = theo_model_ode
% Theophylline oral dose model expressed directly using ODEs

m = InitializeModelAnalytic('TheoModelODEs');

m = AddCompartment(m, 'v', 3, 1);

m = AddSeed(m, 'Xg0', 320); % mg, bolus dose at t = 0 represented as a an initial condition/seed of GI compartment
m = AddSeed(m, 'Xp0', 0); % mg

m = AddState(m, 'Xg', 'v', 'Xg0'); % mg
m = AddState(m, 'Xp', 'v', 'Xp0'); % mg

m = AddParameter(m, 'ka',  2.87); % 1/hr
m = AddParameter(m, 'kel', 0.0723); % 1/hr
m = AddParameter(m, 'V',   32.5); % L

% Swapped order of ODEs to test logic that reorders them to match states
m = AddOde(m, 'Xp', 'ka*Xg - kel*Xp');
m = AddOde(m, 'Xg', '-ka*Xg');

m = AddOutput(m, 'Cp', 'Xp/V'); % mg/L

m = FinalizeModel(m);

end