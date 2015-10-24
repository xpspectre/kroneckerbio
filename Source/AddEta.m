function m = AddEta(m, name, errormodel)
% Add a NONMEM-like inter-individual variability parameter eta to parameter
% name to a Model.Analytic. eta will be called eta_{name}
% Specify error model for name as a function of eta and name. Adds the eta parameter
% and a rule that modifies reaction rates and output so d{name}d{eta} is
% calculated correctly.
%
% Inputs:
%   m: [ Model.Analytic struct ]
%       The model to which the eta will be added
%   name: [ string ]
%       Name of rate parameter that eta provides variability for ("theta")
%   errormodel: [ string {'name*exp(eta)'} ]
%       Expression for error model where an individual's true parameter value is
%       a function of name (theta) and variability (eta)
%
% Outputs:
%   m: [ Model.Analytic struct ]
%       The model with the new eta added.
%
% Note/TODO: currently only allows a single value of eta per parameter. Extend
% to allow multiple for more complicated error models, i.e., additive +
% proportional

% Default eta name
eta = ['eta__' name];

% Generate new ID if not supplied
if nargin < 3
    errormodel = [];
end
if isempty(errormodel)
    errormodel = [name '*exp(' eta ')'];
end

% name and eta can contain invalid chars but aren't quoted because those are names
% Clean up errormodel expression if it contains invalid characters
nameQuoted = name;
if regexp(name, '\W')
    nameQuoted = ['"', nameQuoted, '"'];
end
etaQuoted = eta;
if regexp(eta, '\W')
    etaQuoted = ['"', etaQuoted, '"'];
end
errormodel = replaceSymbolRegex(errormodel, name, nameQuoted);
errormodel = replaceSymbolRegex(errormodel, eta, etaQuoted);

%% Modify parameters and add rule
% Add new parameter for eta
m = AddParameter(m, eta, 0);

% Add rule that replaces k in reaction rates with k w/ error model
m = AddRule(m, [name '__errormodel'], name, errormodel, 'single sub');

m.Ready = false;