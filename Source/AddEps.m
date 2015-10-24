function m = AddEps(m, name, errormodel, eps)
% Add a NONMEM-like intra-individual variability parameter eps to output
% name to a Model.Analytic. Modifies the specified output name to contain the
% error model. A dummy output called Fi_{name} representing the state space output w/o eps is
% added from which to calculate derivatives of eta. eps will be called
% eps__{name}.
%
% Inputs:
%   m: [ Model.Analytic struct ]
%       The model to which the eps will be added
%   name: [ string ]
%       Name of output that eps provides variability for
%   errormodel: [ string {'name + name*eps'} ]
%       Expression for error model where an individual's true output value is
%       a function of name (theta) and variability (eps). Default is
%       proportional model.
%   eps: [ string {eps__{name}} | cell array of strings ]
%       Name of eps matching term in error model. Must specify the cell array of
%       eps names if using a more complicated error model that contains multiple
%       epsilons. Must start with 'eps__'
%
% Outputs:
%   m: [ Model.Analytic struct ]
%       The model with the new eps and modified output with intra-individual
%       error model.

% Standardize inputs
if nargin < 4
    eps = [];
    if nargin < 3
        errormodel = [];
    end
end

% Set default eps name
if isempty(eps) % omitted optional eps name, set to default
    eps = ['eps_' name];
end
% Standardize into cell arrays of string eps names
if ischar(eps)
    eps = {eps};
end
n = length(eps);

% Set default proportional error model
if isempty(errormodel)
    errormodel = [name '*(1 + ' eps{1} ')'];
    
    % Clean up errormodel if eps{1} has an invalid char
    %   Not needed for user supplied errormodels because invalid names should
    %   already be quoted
    if regexp(eps{1}, '\W')
        errormodel = replaceSymbolRegex(errormodel, eps{1}, ['"', eps{1}, '"']);
    end
end

%% Search for output to replace
yNames = {m.Outputs.Name};
yNames(cellfun(@isempty, yNames)) = [];
output = m.Outputs(ismember(yNames, name));
m = RemoveOutput(m, name);

%% Replace output in error model expression with expression for output
% Needed because outputs can't depend on other outputs
% Safe because expression for output is already a valid expression
% Wrap in parentheses to ensure it's added as a unit
errormodel = replaceSymbolRegex(errormodel, name, ['(' output.Expression ')']);

%% Add output for Y (measured solution, with eps, derivatives wrt eps will be taken)
m = AddOutput(m, output.Name, errormodel);

%% Add output for F (state space solution, no eps, derivatives wrt eta will be taken)
m = AddOutput(m, ['Fi__' output.Name], output.Expression);

%% Add new parameters
for i = 1:n
    m = AddParameter(m, eps{i}, 0);
end

m.Ready = false;