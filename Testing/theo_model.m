function m = theo_model(nlme)
% Theophylline oral dose model for nonlinear mixed effects PK/PD modeling
% 1-compartment PK model with 1st order absorption
% Starting param values close to fit values to make tests faster

if nargin < 1
    nlme = [];
end
if isempty(nlme)
    nlme = false;
end

if nlme
    m = InitializeModelAnalytic('TheoModel');
else
    m = InitializeModelAnalytic('TheoModelNlme');
end

m = AddCompartment(m, 'v', 3, 1);

m = AddSeed(m, 'Xg0', 320); % mg, bolus dose at t = 0 represented as a an initial condition/seed of GI compartment
m = AddSeed(m, 'Xp0', 0); % mg

m = AddState(m, 'Xg', 'v', 'Xg0'); % mg
m = AddState(m, 'Xp', 'v', 'Xp0'); % mg

m = AddParameter(m, 'ka',  2.87); % 1/hr (1.5)
m = AddParameter(m, 'kel', 0.0723); % 1/hr (0.15)
m = AddParameter(m, 'V',   32.5); % L (32.5)

m = AddReaction(m, 'Absorption',  {'Xg'}, {'Xp'}, 'ka*Xg'); % don't forget to include reactants in mass-action like analytic rate forms
m = AddReaction(m, 'Elimination', {'Xp'}, {}, 'kel*Xp');

% m = AddOutput(m, 'Xg', 'Xg'); % mg
% m = AddOutput(m, 'Xp', 'Xp'); % mg
m = AddOutput(m, 'Cp', 'Xp/V'); % mg/L

if nlme
    
    % Note: variable names are extracted using regexes which look for parts divided
    % by double-underscores (__), so don't use __ in parameter names above
    
    % Add inter-individual variability parameters etas
    % Remember to add Omega after these
    m = AddEta(m, 'ka'); % default is lognormal, ka = ka*exp(eta_ka), FO method
    m = AddEta(m, 'kel');
    m = AddEta(m, 'V');
    
    % m = AddOmega(m, {'ka','kel','V'});
    m = AddOmega(m, {'ka','kel','V'}, [0.2 0 0; 0.1 0.2 0; 0.1 0.1 0.2]);
    
    % Add intra-individual variability parameters epsilons
    % m = AddErrorModel(m, {'Cp'}, [0.1]); % Simple outer error - constant term only for direct comparison to FO
    m = AddErrorModel(m, {'Cp'}, {'Cp*sigma__Cp1 + sigma__Cp2'}, {'sigma__Cp1', 'sigma__Cp2'}, [0.1, 0.1]); % original combined error
    
end

m = FinalizeModel(m);