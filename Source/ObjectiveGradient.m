function [D, G] = ObjectiveGradient(varargin)
%ObjectiveGradient Evaluate the gradient of a set of objective functions
%
%   This function accepts a FitObject + options or separate model, experimental
%   condition(s), and objective function(s) + options.
%
%   [D, G] = FiniteObjectiveGradient(FitObject, opts)
%   [D, G] = FiniteObjectiveGradient(m, con, obj, opts)
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
%   D [ nT x 1 real vector ]
%       The sum of all objective function gradients
%   G [ scalar real ]
%       The objective function value

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
        assert(nargin >= 3, 'KroneckerBio:ObjectiveGradient:TooFewInputs', 'ObjectiveValue requires at least 3 input arguments')
        m = varargin{1};
        con = varargin{2};
        obj = varargin{3};
        if nargin >= 4
            opts = varargin{4};
        else
            opts = [];
        end
        assert(isscalar(m), 'KroneckerBio:ObjectiveGradient:MoreThanOneModel', 'The model structure must be scalar')
        fit = FitObject.buildFitObject(m, con, obj, opts);
end

%% Run main calculation
[G, D] = fit.computeObjective;
