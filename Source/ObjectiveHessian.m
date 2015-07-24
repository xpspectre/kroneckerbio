function [H, D, G] = ObjectiveHessian(varargin)
%ObjectiveHessian Evaluate the hessian of a set of objective functions
%
%   [H, D, G] = ObjectiveHessian(FitObject, opts)
%   [H, D, G] = ObjectiveHessian(m, con, obj, opts)
%
%   Inputs
%   m: [ model struct scalar ]
%       The KroneckerBio model that will be simulated
%   con: [ experiment struct vector ]
%       The experimental conditions under which the model will be simulated
%   obj: [ objective struct matrix ]
%       The objective structures defining the objective functions to be
%       evaluated.
%       .UseModelSeeds [ logical scalar {false} ]
%           Indicates that the model's seed parameters should be used
%           instead of those of the experimental conditions
%       .UseModelInputs [ logical scalar {false} ]
%           Indicates that the model's inputs should be used instead of
%           those of the experimental conditions
%       .UseParams [ logical vector nk | positive integer vector {1:nk} ]
%           Which kinetic parameters the gradient will be calculated on
%       .UseSeeds [ logical matrix nx by nCon | logical vector nx |
%                   positive integer vector {[]} ]
%           Which seed parameters the gradient will be calculated on
%       .UseControls [ cell vector nCon of logical vectors or positive 
%                      integer vectors | logical vector nq | positive 
%                      integer vector {[]} ]
%           Which input control parameters the gradient will be calculated
%           on
%     	.ObjWeights [ real matrix nObj by nCon {ones(nObj,nCon)} ]
%           Applies a post evaluation weight on each objective function
%           in terms of how much it will contribute to the final objective
%           function value
%       .Normalized [ logical scalar {true} ]
%           Indicates if the gradient should be computed in log parameters
%           space
%    	.UseAdjoint [ logical scalar {false} ]
%           Indicates whether the gradient should be calculated via the
%           adjoint method or the forward method.
%       .RelTol [ nonnegative scalar {1e-6} ]
%           Relative tolerance of the integration
%       .AbsTol [ cell vector of nonnegative vectors | nonnegative vector |
%                 nonegative scalar {1e-9} ]
%           Absolute tolerance of the integration. If a cell vector is
%           provided, a different AbsTol will be used for each experiment.
%       .Verbose [ nonnegative integer scalar {1} ]
%           Bigger number displays more progress information
%
%   Outputs
%       H: [ real vector nT ]
%           The sum of all objective function hessians

% (c) 2013 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

%% Work-up
% Clean up inputs
switch nargin
    case 1
        fit = varargin{1};
    case 2
        fit = varargin{1};
        opts = varargin{2};
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
        assert(isscalar(m), 'KroneckerBio:ObjectiveHessian:MoreThanOneModel', 'The model structure must be scalar')
        fit = FitObject.buildFitObject(m, con, obj, opts);
end

% Put into fit object
fit = FitObject.buildFitObject(m, con, obj, opts);

% For convenience, copy fit object's options into this space
opts = fit.options;

%% Run main calculation
[G, D, H] = fit.computeObjective;

%% Normalization
if opts.Normalized
    T = fit.collectParams;
    nT = length(T);
    
    % Normalize
    H = spdiags(T,0,nT,nT) * H * spdiags(T,0,nT,nT);
    H = full(H); % Matlab bug makes this necessary
    
    D = D .* fit.collectParams;
end