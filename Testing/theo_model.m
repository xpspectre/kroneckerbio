function m = theo_model
% Theophylline oral dose model for nonlinear mixed effects PK/PD modeling
% 1-compartment PK model with 1st order absorption

m = InitializeModelAnalytic('TheoModel');

m = AddCompartment(m, 'v', 3, 1);

m = AddSeed(m, 'Xg0', 320); % mg, bolus dose at t = 0 represented as a an initial condition/seed of GI compartment
m = AddSeed(m, 'Xp0', 0); % mg

m = AddState(m, 'Xg', 'v', 'Xg0'); % mg
m = AddState(m, 'Xp', 'v', 'Xp0'); % mg

m = AddParameter(m, 'ka',  1.5); % 1/hr
m = AddParameter(m, 'kel', 0.15); % 1/hr
m = AddParameter(m, 'V',   25); % L

m = AddReaction(m, 'Absorption',  {'Xg'}, {'Xp'}, 'ka*Xg'); % don't forget to include reactants in mass-action like analytic rate forms
m = AddReaction(m, 'Elimination', {'Xp'}, {}, 'kel*Xp');

% m = AddOutput(m, 'Xg', 'Xg'); % mg
% m = AddOutput(m, 'Xp', 'Xp'); % mg
m = AddOutput(m, 'Cp', 'Xp/V'); % mg/L

% Note: variable names are extracted using regexes which look for parts divided
% by double-underscores (__), so don't use __ in parameter names above
% Add interindividual variability parameters etas
% Remember to add Omega after these
m = AddEta(m, 'ka'); % default is lognormal, ka = ka*exp(eta_ka), FO method
m = AddEta(m, 'kel');
m = AddEta(m, 'V');
m = AddOmega(m, {'ka','kel','V'}, [0.2 0 0; 0.1 0.2 0; 0.1 0.1 0.2]);

% Add intraindividual variability parameters epsilons
% Remember to add Sigma after these
% m = AddEps(m, 'Cp'); % desired default invocation: adds 1 eps for proportional model
m = AddEps(m, 'Cp', 'Cp + Cp*eps__Cp1 + eps__Cp2', {'eps__Cp1', 'eps__Cp2'}); % explicitly provide eps names because of custom error model
m = AddSigma(m, {'Cp1', 'Cp2'}, [0.1 0; 0 0.1]);

m = FinalizeModel(m);

end