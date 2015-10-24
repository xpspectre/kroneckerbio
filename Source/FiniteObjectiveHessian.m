function [H, D, G] = FiniteObjectiveHessian(varargin)
%ObjectiveHessian Evaluate the hessian of a set of objective functions by finite
%   difference approximation
%
%   This function accepts a FitObject + options or separate model, experimental
%   condition(s), and objective function(s) + options.
%
%   [H, D, G] = ObjectiveHessian(FitObject, opts)
%   [H, D, G] = ObjectiveHessian(m, con, obj, opts)
%
%   This function is identical to ObjectiveHessian except that the
%   gradient is calculated differently. This is useful for testing that
%   custom-made objective functions have been coded correctly and that
%   integration tolerances have been set appropriately.
%
% Inputs:
%   FitObject [ scalar FitObject ]
%       Fitting scheme
%   m [ scalar model struct ]
%       The KroneckerBio model that will be simulated
%   con [ nCon x 1 conditions struct vector ]
%       The experimental conditions under which the model will be simulated. A
%       nCon x 1 struct vector is allowed.
%   obj [ nCon x nObj objectives struct matrix ]
%       The objective functions under which the experiments will be simulated.
%       Each row of obj has a corresponding row of con.
%   opts [ options struct ]
%       Options struct allowing the following fields, overriding existing values
%       in m, con, and obj. See FitObject.buildFitObject opts for allowed
%       fields.
%
% Outputs:
%   H [ nT x nT real matrix ]
%       The sum of all objective function hessians
%   D [ nT x 1 real vector ]
%       The sum of all objective function gradients
%   G [ scalar real ]
%       The objective function value

%% Work-up
% Clean up inputs
opts_ = [];
opts_.UseAdjoint = false; % adjoint method can't be used for complex step differentiation for the hessian

switch nargin
    case 1
        fit = varargin{1};
    case 2
        fit = varargin{1};
        opts = varargin{2};
        opts = mergestruct(opts_, opts);
        fit.addOptions(opts);
    otherwise % Old method
        assert(nargin >= 3, 'KroneckerBio:ObjectiveHessian:TooFewInputs', 'ObjectiveValue requires at least 3 input arguments')
        m = varargin{1};
        con = varargin{2};
        obj = varargin{3};
        if nargin >= 4
            opts = varargin{4};
        else
            opts = [];
        end
        opts = mergestruct(opts_, opts);
        assert(isscalar(m), 'KroneckerBio:ObjectiveHessian:MoreThanOneModel', 'The model structure must be scalar')
        fit = FitObject.buildFitObject(m, con, obj, opts);
end

% For convenience, copy fit object's options into this space
opts = fit.options;
verbose = logical(opts.Verbose);

% Store starting parameter sets
T0 = fit.collectParams;

%% Loop through conditions
nT = length(T0);
H = zeros(nT,nT);

% Initial value
if verbose; fprintf('Initial round\n'); end
[G, D] = fit.computeObjective;

for i = 1:nT
    if verbose; fprintf('Step %d of %d\n', i, nT); end
    
    % Set baseline parameters
    T_i = T0(i);
    T_up = T0;
    
    % Change current parameter by finite amount
    step_size = 1e-8;
    if opts.Normalized
        norm_factor = T_i;
    else
        norm_factor = 1;
    end
    if opts.ComplexStep
        imag_factor = 1i;
    else
        imag_factor = 1;
    end
    diff = step_size * norm_factor * imag_factor;
    
    % Compute D
    T_up(i) = T_up(i) + diff;
    fit.updateParams(T_up);
    [~, D_up] = fit.computeObjective;

    % Compute H
    if opts.ComplexStep
        H(:,i) = imag(D_up) ./ step_size;
    else
        H(:,i) = (D_up - D) ./ step_size;
    end
end
