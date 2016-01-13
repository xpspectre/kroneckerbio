function [fitOut, beta, psi, stats] = FitNLME(fit, opts)
% Fit nonlinear mixed effects model using built-in Matlab NLME methods
% Inputs:
%   fit [ FitObject ]
%       Specifies fitting scheme.
%   opts [ options struct ]
%       Options struct allowing the following fields:
%       .Verbose [ scalar nonegative integer {1} ]
%           The amount of displayed information
%       .DisplayTraces [ true | {false} ]
%           Whether to display running plot of parameter estimates while
%           fitting. Useful for seeing convergence in SAEM.
%       .method [ {'LME'} | 'RELME' | 'FO' | 'FOCE' | 'SAEM' ]
%           The estimation/solver method. LME, RELME, FO, and FOCE use nlmefit
%           and SAEM uses nlmefitsa. Warning: FO and FOCE are slow. See the
%           ApproximationType parameter in the nlmefit documentation. Warning:
%           SAEM is slow. Linearized methods are generally faster if they
%           converge.
%       .etas [ cell array of strings ]
%           Names of the model parameters which have inter-individual
%           variability. At least 1 is required.
%       .ErrorModel [ {'constant'} | 'proportional' | 'combined' | 'exponential' ]
%           Form of the error term for the output.
%       .ParamTransform [ 1 x p vector of integers 0-3 {zeros(1,p)} ]
%           Transformation of params, can be integers 0 - 3
%       .FiniteDiffEps [ positive scalar double {1e-6} ]
%           Relative difference used for finite difference gradient. Note: I
%           think original default of eps^(1/3) is too small.
%       .NChains [ positive scalar integer {3} ]
%           Number of Markov chains. SAEM only.
% Outputs:
%   fitOut [ FitObject ]
%       Fit with model updated to fit fixed effects parameters. Conditions and
%       objectives are preserved.
%   beta [ f x 1 double vector ]
%       Estimates for fixed effects parameters
%   psi [ r x r double matrix ]
%       Estimates for the random effects covariate matrix. r = length(etas)
%   stats [ statistics struct ]
%       Usefult stats 
%
% For reference, refer to Mathworks' help pages for nlmefit and nlmefitsa. The
% code in this function uses the same notation.
%   group = subject = patient
%
% Current limitations:
%   Single output/measurement
%   No fitting of inputs or doses; doses specified in the usual way get used as
%       group-specific parameters
% Note: additional parameters for the fitting routines can be specified by
%   directly adding them to the calls below.
% Note: working on sensitivity eq based gradient methods - don't want to use
%   finite differences
% Note: nlmefit uses unconstrained optimization, as does nlmefitsa, so specified
%   bounds for parameters are ignored. log transform params  (ParamTransform = 1) 
%   to ensure they stay positive.
% TODO: expose more options

if nargin < 2
    opts = [];
end

opts_.Verbose         = 1;
opts_.DisplayTraces   = false;
opts_.method          = 'LME';
opts_.etas            = {};
opts_.ErrorModel      = 'constant';
opts_.ParamTransform  = [];
opts_.FiniteDiffEps   = 1e-6;
opts_.NChains         = 3;

opts = mergestruct(opts_, opts);

verbose = opts.Verbose;
if verbose < 0
    error('FitNLME:invalidVerbosity', 'opts.Verbose must be >= 0')
elseif verbose < 1
    statDisplay = 'off';
elseif verbose < 2
    statDisplay = 'final';
else % >= 2
    statDisplay = 'iter';
end

method = opts.method;
if ~(strcmp(method, 'LME') || strcmp(method, 'RELME') || strcmp(method, 'FO') || strcmp(method, 'FOCE') || strcmp(method, 'SAEM'))
    error('FitNLME:invalidMethod', 'Method %s not recognized', method)
end

errorModel = opts.ErrorModel;
if ~(strcmp(errorModel, 'constant') || strcmp(errorModel, 'proportional') || strcmp(errorModel, 'combined') || strcmp(errorModel, 'exponential'))
    error('FitNLME:invalidErrorModel', 'Error model %s not recognized', errorModel)
end

%% Extract a common model and condition for which calculations will be performed on
% Assumes that all groups/subjects are variants of this single model/condition
% Update* is called on these to efficiently calculate values
model = fit.Models(fit.componentMap(1,1));
condition = fit.Conditions(fit.componentMap(1,2));

%% Get output index
output = ismember({model.Outputs.Name}, fit.Objectives(1).Outputs);

%% Get parameters
etas = opts.etas;
assert(iscell(etas), 'FitNLME:invalidEtas', 'etas is not a cell array')

% Extract fit parameter names and indices from model/condition
k = model.k;
kFit = logical(condition.ParamsSpec{1})';
kNames = {model.Parameters.Name};
kNames = kNames(kFit);

% Fixed effects/all fit parameters
p = sum(kFit); % number of parameters
f = p; % number of fixed effects
fixedEffects = true(1,f); % assume all parameters have fixed effects

% Random effects/inter-individual variability
randomEffects = ismember(kNames, etas); % only specified parameters have random effects
assert(sum(randomEffects) == length(etas), 'FitNLME:numEtasMismatch', 'Number of found etas %g not equal to number of specified etas %g', sum(randomEffects), length(etas))

%% Extract measurements from fit
% Each group corresponds to the measurements from a single objective function
%   Same order as specified in fit
m = fit.nObjectives; % number of groups (subjects)
times        = [];
measurements = [];
group        = [];
ni = zeros(1,m); % number of observations in group i
for i = 1:m
    times_i        = vec(fit.Objectives(i).Times);
    measurements_i = vec(fit.Objectives(i).Measurements);
    ni_i           = length(times_i);
    assert(ni_i == length(measurements_i), 'FitNLME:measurementsMismatch', 'Number of times and measurements don''t match for subject %g', i)
    group_i        = i*ones(ni_i, 1);
    
    times        = [times; times_i];
    measurements = [measurements; measurements_i];
    group        = [group; group_i];
    ni(i)        = ni_i;
end
n = length(times); % number of observations
h = 1; % single predictor variable time; put inputs, doses, etc. in group-specific variables V

% Hack to convert 0's and small negatives into +eps. Needed for exponential
% error model
measurements(measurements<=0) = 1e-3;

%% Group-specific predictors (subject-specific, covariates)
% Includes kroneckerbio s, q, and h
% All scalars
% All subjects must have same specified group-specific predictors, i.e., all must have a value for a dose
% Vmap relates V to s, q, and h, where s = 1, q = 2, and h = 3
V = [];
for i = 1:m
    con = fit.Conditions(i);
    V = [V; con.s', con.q', con.h'];
    if i == 1 % only need one, as in assumption above
        Vmap = [ones(1,con.ns), 2*ones(1,con.nq), 3*ones(1,con.nh)];
    end
end
g = length(Vmap); % number of group-specific predictor variables

%% Function handle for model predictions
fun = @modelfun;

%% Extract starting parameter values
beta0 = model.k(kFit);
% beta0 = [1 1 1]'; % test new starting point
assert(length(beta0) == p, 'FitNLME:numBetasMismatch', 'Number of found initial guesses for fixed effects parameters %g not equal to number of specified fixed effects parameters %g', length(beta0), p)

% Standardize paramTransform
paramTransform = opts.ParamTransform;
if isempty(paramTransform)
    paramTransform = zeros(size(beta0))';
elseif isscalar(paramTransform)
    paramTransform = paramTransform*ones(size(beta0))';
end

%% Additional options
statOpts = statset('Display', statDisplay);
if opts.DisplayTraces
    statOpts.OutputFcn = @nlmefitoutputfcn;
end

%% Fit using specified method
if strcmp(method, 'LME') || strcmp(method, 'RELME') || strcmp(method, 'FO') || strcmp(method, 'FOCE')
    [beta, psi, stats] = nlmefit(times, measurements, group, V, fun, beta0, ...
        'FEParamsSelect', fixedEffects, ...
        'REParamsSelect', randomEffects, ...
        'Vectorization', 'SingleGroup', ...  % each subject is integrated once for all times
        'ErrorModel', errorModel, ...
        'ParamTransform', paramTransform, ...
        'Options', statOpts, ...
        'OptimFun', 'fminunc', ...
        'ApproximationType', method);
elseif strcmp(method, 'SAEM')
    [beta, psi, stats] = nlmefitsa(times, measurements, group, V, fun, beta0, ...
        'FEParamsSelect', fixedEffects, ...
        'REParamsSelect', randomEffects, ...
        'Vectorization', 'SingleGroup', ...  % each subject is integrated once for all times
        'ErrorModel', errorModel, ...
        'ParamTransform', paramTransform, ...
        'Options', statOpts, ...
        'OptimFun', 'fminunc', ...
        'NBurnIn', 3, ...
        'NChains', opts.NChains);
end

%% Make FitObject with fit results
% Contains individual conditions and objectives but a single model with averaged
% values.
fitOut = fit;
for i = 1:fitOut.nModels
    kNew = k;
    kNew(kFit) = beta;
    
    m = fit.Models(i).Update(kNew);
    m.BaseModel = fit.Models(i).BaseModel; % need this to have similar structs
    
    fitOut.Models(i) = m;
end

%% Helper functions
    function yfit = modelfun(phi, xfun, vfun)
        % Model function that returns predicted values to compare to
        % measurements.
        % Inputs:
        %   phi: 1 x p vector of model parameters
        %   xfun: l x 1 vector of times for just this group
        %   vfun: 1 x g vector of group-specific predictors for a single group
        % Outputs:
        %   yfit: l x 1 vector of predicted values, 1 for each xfun/time
        %
        % Sanity checks with asserts
        % Note: need to be able to exit and retry w/ new values
        assert(all(size(phi) == [1,p]), 'FitNLME:modelfun:phiDims', 'Invalid dims for phi')
        assert(size(xfun,2) == 1, 'FitNLME:modelfun:xfunDims', 'Invalid dims for xfun')
        assert(all(size(vfun) == [1,g]), 'FitNLME:modelfun:vfunDims', 'Invalid dims for vfun')
        
        % Specify final time for integration
        %   xfun -> times
        tf = xfun(end);
        
        % Set up model and condition for given arguments
        %   phi -> k in m
        %   vfun -> s, q, and h in con
        kNew = k;
        kNew(kFit) = phi;
        sNew = vfun(Vmap==1);
        qNew = vfun(Vmap==2);
        hNew = vfun(Vmap==3);
        
        model = model.Update(kNew);
        condition = condition.Update(sNew, qNew, hNew);
        
        intOpts = [];
        intOpts.Verbose = 1; % 0 for no output, 1 to display a line for each integration
        sim = SimulateSystem(model, condition, tf, intOpts);
        yfit = sim.y(xfun, output)';
        
        assert(length(yfit) == length(xfun), 'FitNLME:modelfun:yfitDims', 'yfit dims don''t match xfun dims')
    end

end