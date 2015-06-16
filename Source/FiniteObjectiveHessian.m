function [H, All] = FiniteObjectiveHessian(m, con, obj, opts)
%FiniteObjectiveHessian Compute the hessian of an objective function by
%   the finite difference approximation
%
%   Mathematically: H = dG/dp2 or H = dG/dlnp2
%
%   [H, All] = FiniteObjectiveHessian(m, con, obj, opts)
%
%   Inputs:
%       m    - The KroneckerBio model that will be simulated
%       con  - A structure vector of the experimental conditions under
%              which the hessian will be evaluated
%       obj  - A structure array of the objective functions under which the
%              hessian will be evaluated. Each row of obj is matched to the
%              corresponding entry in con
%       opts - Optional function options
%           .useParams   - Vector of indexes indicating the rate constants
%                          whose sensitivities will be considered
%           .useICs      - Vector of indexes indicating the initial
%                          concentrations whose sensitivities will be
%                          considered
%           .UseModelICs - Boolean that determines whether to use the
%                          initial concentrations of the model or the
%                          conditions. Default = true
%           .Normalized  - Boolean that determines if the simple hessian or
%                          the normalized hessian will be computed. The
%                          normalized hessian is normalized with respect to
%                          the values of the parameters. Default = true
%           .Verbose     - Print progress to command window
%   Outputs:
%       H = FiniteObjectiveHessian(m, con, obj, ...)
%           H - The hessian as an array
%
%       [H, All] = FiniteObjectiveHessian(m, con, obj, ...)
%           All - The hessian is the sum of hessians from all experiements,
%           but All returns the hessians in a size(obj) cell array
%
%   Additional info:
%   - The experimental condition vector can also be a cell vector
%   - The objective function array can also be a cell array. Empty entries 
%   in the cell array and entries in the structure array with empty
%   values are ignored. This way, conditions can have different numbers of
%   objective functions associated with them.

% (c) 2015 David R Hagen and Bruce Tidor
% This software is released under the GNU GPLv2 or later.

%% Work-up
% Clean up inputs
if nargin < 4
    opts = [];
end

assert(nargin >= 3, 'KroneckerBio:ObjectiveGradient:TooFewInputs', 'ObjectiveGradient requires at least 3 input arguments')
assert(isscalar(m), 'KroneckerBio:ObjectiveGradient:MoreThanOneModel', 'The model structure must be scalar')

% Put into fit object
fit = FitObject.buildFitObject(m, con, obj, opts);

% For convenience, copy fit object's options into this space
opts = fit.options;
verbose = logical(opts.Verbose);

% Store starting parameter sets
T0 = fit.collectParams;

%% Loop through conditions
nT = length(T0);
H = zeros(nT,1);

% Initial value
if verbose; fprintf('Initial round\n'); end
[~, D] = fit.computeObjective;

for i = 1:nT
    if verbose; fprintf('Step %d of %d\n', i, nT); end
    
    % Set baseline parameters
    T_i = T0(i);
    T_up = T0;
    
    % Change current parameter by finite amount
    if opts.Normalized
        diff = T_i * 1e-8;
    else
        diff = 1e-8;
    end
    
    % Compute objective values
    T_up(i) = T_up(i) + diff;
    fit.updateParams(T_up);
    [~, D_up] = fit.computeObjective;

    % Compute D
    if opts.Normalized
        H(:,i) = T_i * T0 .* (D_up - D) ./ diff;
    else
        H(:,i) = (D_up - D) ./ diff ;
    end
end
