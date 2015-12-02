function m = AddErrorModel(m, name, errormodel, sigma, guess)
% Add intra-individual variability / error model R to Model.Analytic for output.
% The error is eps ~ N(0,R) distributed.
%   The output is renamed h__{name}
%   A dummy output Ri__{name} is added to denote the error model
%
% Inputs:
%   m: [ Model.Analytic struct ]
%       The model to which the error model will be added
%   name: [ string ]
%       Name of output that eps provides variability for
%   errormodel: [ string {'sigma__{name}'} ]
%       Expression for error model where an individual's true output value is
%       a function of name (theta) and variability (sigma). Default is a
%       constant error model.
%   sigma: [ string {sigma__{name}} | cell array of strings ]
%       Name of sigma(s) matching term in error model. Must specify the cell array of
%       sigma names if using a more complicated error model that contains multiple
%       sigmas like a combined constant + proportional model. Must start with
%       'sigma__'.
%   guess: [ nsigma array of dobules {0.1} ]
%       Initial guesses for sigma, matching number and order of sigmas
%
% Outputs:
%   m: [ Model.Analytic struct ]
%       The model with the new sigma and modified outputs

% Standardize inputs
if nargin < 5
    guess = [];
    if nargin < 4
        sigma = [];
        if nargin < 3
            errormodel = [];
        end
    end
end

% Set default eps name
if isempty(sigma) % omitted optional sigma name, set to default
    sigma = ['sigma__' name];
end

% Standardize into cell arrays of string sigma names
if ischar(sigma)
    sigma = {sigma};
end
n = length(sigma);

% Standardize initial guesses
if isempty(guess)
    guess = repmat(0.1, [n,1]);
elseif numel(guess) ~= n
    error('AddErrorModel:WrongNumberOfGuesses', 'Number of guesses doesn''t match number of sigmas')
end

% Set default constant error model
if isempty(errormodel)
    errormodel = sigma{1};
    
    % Clean up errormodel if sigma{1} has an invalid char
    %   Not needed for user supplied errormodels because invalid names should
    %   already be quoted
    if regexp(sigma{1}, '\W')
        errormodel = replaceSymbolRegex(errormodel, sigma{1}, ['"', sigma{1}, '"']);
    end
end

%% Search for output to replace
yNames = {m.Outputs.Name};
yNames(cellfun(@isempty, yNames)) = [];

y = ismember(yNames, name);
if any(y) % Pull output existing in finalized model
    output = m.Outputs(y);
    m.Outputs(y) = [];
    m.ny = m.ny - 1;
else % Add completely new output
    error('AddErrorModel:InvalidOutput', 'Output %s for intra-individual variability not found in model.', name)
end

%% Replace output in error model expression with expression for output
% Needed because outputs can't depend on other outputs
% Safe because expression for output is already a valid expression
% Wrap in parentheses to ensure it's added as a unit
% Square the entire term to fit Ri form
errormodel = ['(' replaceSymbolRegex(errormodel, name, ['(' output.Expression ')']) ')^2'];

%% Add output for Y (measured solution, with eps, derivatives wrt eps will be taken)
m = AddOutput(m, ['Ri__' output.Name], errormodel);

%% Add output for h (state space solution, no eps, derivatives wrt eta will be taken)
m = AddOutput(m, ['h__' output.Name], output.Expression);

%% Add new parameters
for i = 1:n
    m = AddParameter(m, sigma{i}, guess(i)); % some random value
end

m.Ready = false;