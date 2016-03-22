function [F, All] = ObjectiveInformation(varargin)
%ObjectiveInformation Compute the Fisher information matrix (FIM) of a set of
%   objective functions
%
%   Mathematically: F = E(d2logpdT2(y,Y), y, p_y|T(y,T))
%
%   [F, All] = ObjectiveInformation(FitObject)
%   [...] = ObjectiveInformation(FitObject, opts)
%   [...] = ObjectiveInformation(FitObject, opts, dxdTSol)
%   [...] = ObjectiveInformation(m, con, obj)
%   [...] = ObjectiveInformation(m, con, obj, opts)
%   [...] = ObjectiveInformation(m, con, obj, opts, dxdTSol)
%
% Inputs:
%   FitObject [ scalar FitObject ]
%       Fitting scheme
%   m [ scalar model struct ]
%       The KroneckerBio model that will be simulated
%   con [ nCon x 1 conditions struct vector ]
%       The experimental conditions under which the model will be simulated
%   obj [ nCon x nObj objectives struct matrix ]
%       The objective structures defining the objective functions to be
%       evaluated. Each row of obj has a corresponding row of con.
%   opts [ options struct ]
%       Options struct allowing the following fields, overriding existing values
%       in m, con, and obj. See FitObject.buildFitObject opts for allowed
%       fields.
%   dxdTSol [ nCon x 1 cell array of sensitivity solution structs | sensitivity solution struct {} ]
%       A structure vector containing the solution to the model
%       sensitivities under each condition can be provided to prevent this
%       method from recalculating it
%
% Outputs:
%   F [ double matrix ]
%       Combined fisher information matrix. The FIM is the sum of all
%       fisher information matrices assuming there is no
%       covariance between errors in seperate experiments.
%   All [ nCon x 1 cell array of double matrices ]
%       The individual FIMs for each experiment

%% Work-up
% Clean up inputs
switch nargin
    case 1
        fit = varargin{1};
    case 2
        fit = varargin{1};
        opts = varargin{2};
        fit.addOptions(opts);
    otherwise
        if isa(varargin{1}, 'FitObject')
            fit = varargin{1};
            opts = varargin{2};
            fit.addOptions(opts);
            dxdTSol = varargin{3};
        else
            assert(nargin >= 3, 'KroneckerBio:ObjectiveInformation:TooFewInputs', 'ObjectiveInformation requires at least 3 input arguments')
            m = varargin{1};
            con = varargin{2};
            obj = varargin{3};
            if nargin >= 4
                opts = varargin{4};
            else
                opts = [];
            end
            if nargin >= 5
                dxdTSol = varargin{5};
            else
                dxdTSol = [];
            end
            assert(isscalar(m), 'KroneckerBio:ObjectiveInformation:MoreThanOneModel', 'The model structure must be scalar')
            fit = FitObject.buildFitObject(m, con, obj, opts);
        end
end

nCon = fit.nConditions;

% Clean up dxdTSol
if exist('dxdTSol', 'var') && ~isempty(dxdTSol)
    if ~iscell(dxdTSol)
        dxdTSol = {dxdTSol};
    else
        assert(all(size(dxdTSol) == [nCon,1]), 'KroneckerBio:ObjectiveInformation:InvalidSensitivity', 'dxdTSol must match the supplied conditions.')
    end
else
    dxdTSol = {};
end

%% Run main calculation
[F, All] = fit.computeInformation(dxdTSol);
