function m = AddErrorModel(varargin)
% Add intra-individual variability / error model R to Model.Analytic for output.
% The error is eps ~ N(0,R) distributed.
%   The output is renamed h__{name}
%   A dummy output Ri__{name} is added to denote the error model
%   Additional 
%
% Usage:
%   m = AddErrorModel(m, names)
%   m = AddErrorModel(m, names, errorvalues)
%   m = AddErrorModel(m, names, errormodel, errornames)
%   m = AddErrorModel(m, names, errormodel, errornames, errorvalues)
%
% Inputs:
%   m: [ Model.Analytic struct ]
%       The model to which the error model will be added
%   names: [ cell array of strings ]
%       Name of outputs that with measurements that error model applies to
%   errormodel: [ cell matrix of strings {'sigma__{name}'} ]
%       Error model expressions. Symbolically evaluated to a covariance matrix.
%       Default is constant, independent error for each output.
%   errornames: [ nSigma cell array of strings {sigma__{name} ]
%       Name of sigma(s) matching term in error model. Must specify the cell array of
%       sigma names if using a more complicated error model that contains multiple
%       sigmas like a combined constant + proportional model. Must start with
%       'sigma__'.
%   errorvalues: [ nSigma double array {} ]
%       Initial guesses for sigmas, matching number and order of sigmas. Note:
%       can't set default values to 10% of output like for Omega because we need
%       the symbolic expression first. TODO: consider adding this functionality
%       to FinalizeModel
%
% Outputs:
%   m: [ Model.Analytic struct ]
%       The model with the new sigma and modified outputs

% Standardize inputs
switch nargin
    case 2
        m           = varargin{1};
        names       = varargin{2};
        errormodel  = [];
        errornames  = [];
        errorvalues = [];
    case 3
        m           = varargin{1};
        names       = varargin{2};
        errormodel  = [];
        errornames  = [];
        errorvalues = varargin{3};
    case 4
        m           = varargin{1};
        names       = varargin{2};
        errormodel  = varargin{3};
        errornames  = varargin{4};
        errorvalues = [];
    case 5
        m           = varargin{1};
        names       = varargin{2};
        errormodel  = varargin{3};
        errornames  = varargin{4};
        errorvalues = varargin{5};
    otherwise
        error('AddErrorModel:InvalidNumberArgs', 'Invalid number of arguments (%g) specified', nargin)
end

names = vec(names);
nOutputs = length(names);

% Set default constant error model if not specified
if isempty(errormodel) % implies errornames is also empty, but can still have initial guesses assuming default constant error model
    errormodel = cell(nOutputs);
    errornames = cell(nOutputs,1);
    for i = 1:nOutputs
        errorname = ['sigma__' names{i}];
        errornames{i} = errorname;
        
        % Quote errorname in errormodel if name has an invalid char
        if regexp(errorname, '\W')
            errorname = ['"', errorname, '"'];
        end
        errormodel{i,i} = errorname;
    end
end
nSigma = length(errornames);

% Verify errormodel is the right shape
assert(all(size(errormodel) == [nOutputs,nOutputs]), 'AddErrorModel:InvalidErrorModelShape', 'errormodel must be a nName x nName cell matrix. Specified error model was %g x %g', size(errormodel,1), size(errormodel,2))

% TODO: verify errormodel is lower triangular

% Verify that error models refer to outputs that exist in model
y_names = vec({m.Outputs.Name});
y_names = y_names(~cellfun('isempty', y_names));
foundNameMask = ismember(names, y_names);
assert(all(foundNameMask), 'AddErrorModel:InvalidName', 'Names: %s not found in model', cellstr2str(names(~foundNameMask)))

% Verify errornames is a cell array of strings
assert(iscellstr(errornames), 'AddErrorModel:InvalidErrorNames', 'errornames must be a cell array of strings')

% Set default values of sigmas if not given
if isempty(errorvalues)
    warning('AddErrorModel:ErrorValuesNotSpecified', 'Sigma values not specified. Setting initial guesses to 0.1. Note/TODO: implement more sane guess after model is finalized')
    errorvalues = repmat(0.1, nSigma, 1);
end

% Verify number of guesses matches number of names specified
% TODO: verify that errornames and errorvalues match errormodel - more difficult
% - requires symbolic substitutions
assert(numel(errornames) == numel(errorvalues), 'AddErrorModel:WrongNumberOfErrorValues', 'Number of error values doesn''t match number of error names')

%% Create new outputs for errormodels
% Needed because outputs can't depend on other outputs
% Safe because expression for output is already a valid expression
% Wrap in parentheses to ensure it's added as a unit
% Square the entire term to fit Ri form
% Only operate on lower diagonal
% Only add nonempty error models (i.e., errormodel terms with parameters, ignoring covariances = 0 usually) because
%   calculating extra sensitivities scales badly - TODO: refactor sensitivity eq
%   assembly to not have to do these extra terms
y_exprs = vec({m.Outputs.Expression});
y_exprs = y_exprs(~cellfun('isempty', y_exprs));

for i = 1:nOutputs
    for j = 1:i
        em = errormodel{i,j};
        if ~isempty(em)
            for k = 1:nOutputs
                em = replaceSymbolRegex(em, y_names{k}, ['(' y_exprs{k} ')']);
            end
            em = ['(' em ')^2'];
            m = AddOutput(m, ['Ri__' y_names{i} '__' y_names{j}], em);
        end
    end
end

%% Add new parameters
for i = 1:nSigma
    m = AddParameter(m, errornames{i}, errorvalues(i));
end

m.Ready = false;