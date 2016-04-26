function [m, con, obj, options] = ConvertFitOpts(m, con, obj, opts)
% Convert FitObjective's old input args into new input args format
%   opts: [ options struct scalar {} ]
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
%           Indicates the input control parameters that will be allowed to
%           vary during the optimization
%       .UseDoseControls [ cell vector nCon of logical vectors or positive 
%                           integer vectors | logical vector nq | positive 
%                           integer vector {[]} ]
%           Indicates the input control parameters that will be allowed to
%           vary during the optimization
%       .LowerBound [ nonegative vector {0} ]
%           The lower bound on the fitted parameters. It can be length
%           nk+n_con*nx, nk+nx, nT, just nk if nTx = 0, or a scalar. The
%           bounds will be interpreted in that order if the length matches
%           multiple orders.
%       .UpperBound [ nonegative vector {0} ]
%           The upper bound for the fitted parameters. It must be the same
%           length as LowerBound.
%     	.ObjWeights [ real matrix nCon by nObj {ones(nCon,nObj)} ]
%           Applies a post evaluation weight on each objective function
%           in terms of how much it will contribute to the final objective
%           function value.
%       .Normalized [ logical scalar {true} ]
%           Indicates if the optimization should be done in log parameters
%           space
%    	.UseAdjoint [ logical scalar {false} ]
%           Indicates whether the gradient should be calculated via the
%           adjoint method or the forward method
%     	.TolOptim [ positive scalar {1e-5} ]
%           The objective tolerance. The optimization stops when it is
%           predicted that the objective function cannot be improved more
%           than this in the next iteration.
%     	.Restart [ nonnegative integer scalar {0} ]
%           A scalar integer determining how many times the optimzation
%           should restart once optimization has stopped.
%     	.RestartJump [ handle @(iter,G) returns nonnegative vector nT or
%                      scalar | nonnegative vector nT or scalar {0.001} ]
%           This function handle controls the schedule for the noise that
%           will be added to the parameters before each restart. The
%           parameters for the next iteration will be normally distributed
%           in log space with a mean equal to the previous iteration and a
%           standard deviation equal to the value returned by this
%           function. The value returned should usually be a scalar, but it
%           can also be a vector with length equal to the number of active
%           parameters. It can also be numeric, and the noise will be
%           treated as this constant value.
%      	.TerminalObj [ real scalar {-inf} ]
%           Optimization is halted when this objective function value is
%           reached
%       .MaxStepSize [ nonegative scalar {1} ]
%           Scalar fraction indicator of the maximum relative step size
%           that any parameter can take in a single interation
%     	.Algorithm [ string {active-set} ]
%           Option for fmincon. Which optimization algorithm to use
%     	.MaxIter [ postive scalar integer {1000} ]
%           Option for fmincon. Maximum number of iterations allowed before
%           optimization will be terminated.
%     	.MaxFunEvals [ postive scalar integer {5000} ]
%           Option for fmincon. Maximum number of objective function
%           evaluations allowed before optimization will be terminated.
%       .RelTol [ nonnegative scalar {1e-6} ]
%           Relative tolerance of the integration
%       .AbsTol [ nCon x 1 cell vector of nonnegative vectors | nonnegative vector |
%                 nonegative scalar {1e-9} ]
%           Absolute tolerance of the integration. If a cell vector is
%           provided, a different AbsTol will be used for each experiment.
%       .Verbose [ nonnegative integer scalar {1} ]
%           Bigger number displays more progress information
%       .GlobalOptimization [ logical scalar {false} ]
%           Use global optimization in addition to fmincon
%       .GlobalOpts [ options struct scalar {} ]
%           TODO: API in progress
%
% Note:
% - The options above may not be complete. Some new functionality may be
% allowed. See BuildFitOpts for all allowed options.

options = BuildFitOpts;

if nargin == 4
    options = BuildFitOpts(options, opts);
end

assert(isscalar(m), 'KroneckerBio:ConvertFitOpts:MoreThanOneModel', 'The model structure must be scalar') % enforce old input args restriction

% Ensure structures are proper sizes
[con, nCon] = fixCondition(con);
[obj, nObj] = fixObjective(obj, nCon);

nk = m.nk;
ns = con(1).ns; % required all conditions have same number of seeds - could easily remove this restriction
nqs = [con.nq];
nhs = [con.nh];

% Ensure UseParams is logical vector
if ~isfield(opts, 'UseParams')
    opts.UseParams = true(nk,1);
end
[opts.UseParams, nTk] = fixUseParams(opts.UseParams, nk);

% Ensure UseSeeds is a logical matrix
if ~isfield(opts, 'UseSeeds')
    opts.UseSeeds = false(ns,1);
end
[opts.UseSeeds, nTs] = fixUseSeeds(opts.UseSeeds, ns, nCon);

% Ensure UseControls is a cell vector of logical vectors
if ~isfield(opts, 'UseInputControls')
    for iCon = 1:nCon
        opts.UseInputControls = false(nqs(iCon),1);
    end
end
if ~isfield(opts, 'UseDoseControls')
    for iCon = 1:nCon
        opts.UseDoseControls = false(nhs(iCon),1);
    end
end
[opts.UseInputControls, nTq] = fixUseControls(opts.UseInputControls, nCon, cat(1,con.nq));
[opts.UseDoseControls, nTh] = fixUseControls(opts.UseDoseControls, nCon, cat(1,con.nh));

if isfield(opts, 'ObjWeights')
    assert(all(size(opts.ObjWeights) == [nCon,nObj]))
else
    opts.ObjWeights = ones(nCon,nObj);
end

[opts.continuous, opts.complex, opts.tGet] = fixIntegrationType(con, obj);

for iCon = 1:nCon
    opts_i = [];
    
    % Convert UseParams to expected form; old form is nk x 1 logical vector
    usek = zeros(nk,1);
    usek(opts.UseParams) = 1:sum(opts.UseParams);
    opts_i.UseParams = usek;
    
    % Convert UseSeeds to expected form; old form is ns x nCon logical matrix
    uses = zeros(ns,1);
    uses(opts.UseSeeds(:,iCon)) = 1:sum(opts.UseSeeds(:,iCon));
    opts_i.UseSeeds = uses;
    
    % Convert UseInputControls to expected form
    useq = zeros(nqs(iCon),1);
    useq(opts.UseInputControls{iCon}) = 1:sum(opts.UseInputControls{iCon});
    opts_i.UseInputControls = useq;
    
    % Convert UseDoseControls to expected form
    useh = zeros(nhs(iCon),1);
    useh(opts.UseDoseControls{iCon}) = 1:sum(opts.UseDoseControls{iCon});
    opts_i.UseDoseControls = useh;
    
    opts_i.ObjWeights = opts.ObjWeights(iCon,:)';
    opts_i.Continuous = opts.continuous(iCon);
    opts_i.Complex = opts.complex(iCon);
    opts_i.tGet = opts.tGet{iCon};
    
    options = BuildFitOpts(options, m, con(iCon), obj(iCon,:), opts_i);
end

obj = vec(obj); % rearrangned to be the right format; m and con are already good