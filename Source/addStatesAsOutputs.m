function m = addStatesAsOutputs(m, UseShortNames)
%addStatesAsOutputs A quick and dirty helper script that adds one output
%   for each state currently in the model. Note: this
%   function is order-dependent, meaning it only matches species already in the
%   model.
%
%   m = addStatesAsOutputs(m)
%
% Inputs:
%   m [ model ]
%       The model whose existing states are added as outputs
%   UseShortNames [ true | {false} ]
%       Whether to use just the species names without compartment prefixes.
%       Throws an error if multiple species have the same name (for multiple
%       compartments)
%
% Outputs:
%   m [ model ]
%       Model with states added as outputs

if nargin < 2
    UseShortNames = false;
end

assert(isscalar(UseShortNames) && islogical(UseShortNames), 'KroneckerBio:addStatesAsOutputs:InvalidUseShortNamesOption', 'useShortNames arg must be a scalar logical')

% Note that using short names for the expressions should be safe if they're
%   already validated to be unique in the model
if UseShortNames
    names = vec({m.States(1:m.nx).Name});
    assert(numel(names) == numel(unique(names)), 'KroneckerBio:addStatesAsOutputs:NonUniqueSpecies', 'useShortNames was set to true but multiple species share the same name')
else
    names = vec(strcat({m.States(1:m.nx).Compartment}, '.', {m.States(1:m.nx).Name}));
end

for i = 1:numel(names)
    if is(m, 'Model.MassActionAmount')
        m = AddOutput(m, names{i});
    elseif is(m, 'Model.Analytic')
        m = AddOutput(m, names{i}, ['"' names{i} '"']); % quotes around expressions with potentially invalid names
    else
        error('KroneckerBio:AddState:m', 'm must be a model')
    end
end
