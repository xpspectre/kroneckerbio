function m = theo_model_builtin
% Model for built-in Matlab fitnlme and fitnlmesa functions. Doesn't require
% augmenting model with error and variance terms
% Theophylline oral dose model for nonlinear mixed effects PK/PD modeling
% 1-compartment PK model with 1st order absorption

m = InitializeModelAnalytic('TheoModel');

m = AddCompartment(m, 'vol', 3, 1);

m = AddSeed(m, 'Xg0', 320); % mg, bolus dose at t = 0 represented as a an initial condition/seed of GI compartment
m = AddSeed(m, 'Xp0', 0); % mg

m = AddState(m, 'Xg', 'vol', 'Xg0'); % mg
m = AddState(m, 'Xp', 'vol', 'Xp0'); % mg

m = AddParameter(m, 'ka',  1.5); % 1/hr
m = AddParameter(m, 'kel', 0.15); % 1/hr
m = AddParameter(m, 'V',   25); % L

m = AddReaction(m, 'Absorption',  'Xg', '', 'Xp', '', 'ka*Xg'); % don't forget to include reactants in mass-action like analytic rate forms
m = AddReaction(m, 'Elimination', 'Xp', '', '',   '', 'kel*Xp');

% m = AddOutput(m, 'Xg', 'Xg'); % mg
% m = AddOutput(m, 'Xp', 'Xp'); % mg
m = AddOutput(m, 'Cp', 'Xp/V'); % mg/L

m = FinalizeModel(m);

end