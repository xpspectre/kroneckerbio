function G = ObjectiveValue(varargin)
%ObjectiveValue Evaluate a set of objective functions
%
%   G = ObjectiveValue(FitObject, opts)
%   G = ObjectiveValue(m, con, obj, opts)
%   
%   Inputs
%   m: [ model struct scalar ]
%       The KroneckerBio model that will be simulated
%   con: [ experiment struct vector ]
%       The experimental conditions under which the model will be simulated
%   obj: [ objective struct matrix ]
%       The objective structures defining the objective functions to be
%       evaluated.
%   opts: [ options struct scalar ]
%       Optional
%       .UseParams [ logical vector nk | positive integer vector {1:nk} ]
%           Indicates the kinetic parameters that will be allowed to vary
%           during the optimization
%       .UseSeeds [ logical matrix ns by nCon | logical vector ns |
%                   positive integer vector {[]} ]
%           Indicates the seeds that will be allowed to vary during the
%           optimzation. If UseModelSeeds is true then UseSeeds can be a
%           vector of linear indexes or a vector of logicals length of ns.
%           If UseModelSeeds is false then UseSeeds can be a matrix of
%           logicals size ns by nCon. It can also be a vector of length ns,
%           and every experiment will be considered to have the same active
%           seed parameters. It can also be a vector of linear indexes into
%           the ns vector and assumed the same for all conditions.
%       .UseInputControls [ cell vector nCon of logical vectors or positive 
%                           integer vectors | logical vector nq | positive 
%                           integer vector {[]} ]
%           Indicates the input control parameters that are active
%       .UseDoseControls [ cell vector nCon of logical vectors or positive 
%                           integer vectors | logical vector nq | positive 
%                           integer vector {[]} ]
%           Indicates the dose control parameters that are active
%       .RelTol [ nonnegative scalar {1e-6} ]
%           Relative tolerance of the integration
%       .AbsTol [ cell vector of nonnegative vectors | nonnegative vector |
%                 nonegative scalar {1e-9} ]
%           Absolute tolerance of the integration. If a cell vector is
%           provided, a different AbsTol will be used for each experiment.
%     	.ObjWeights [ real matrix nObj by nCon {ones(nObj,nCon)} ]
%           Applies a post evaluation weight on each objective function
%           in terms of how much it will contribute to the final objective
%           function value.
%       .Verbose [ nonnegative integer scalar {1} ]
%           Bigger number displays more progress information
%
%   Outputs
%       G: [ real scalar ]
%           The sum of all objective function values

% (c) 2015 David R Hagen & Bruce Tidor
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
        assert(nargin >= 3, 'KroneckerBio:ObjectiveValue:TooFewInputs', 'ObjectiveValue requires at least 3 input arguments')
        m = varargin{1};
        con = varargin{2};
        obj = varargin{3};
        if nargin >= 4
            opts = varargin{4};
        else
            opts = [];
        end
        assert(isscalar(m), 'KroneckerBio:ObjectiveValue:MoreThanOneModel', 'The model structure must be scalar')
        fit = FitObject.buildFitObject(m, con, obj, opts);
end

%% Run appropriate objective evaluation
G = fit.computeObjective;
