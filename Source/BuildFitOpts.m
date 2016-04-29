function opts = BuildFitOpts(opts, varargin)
% Build options for advanced fitting schemes. Each set of
%   model+condition+objective corresponds to an experimental condition, it's
%   dependent model, and all the objective functions for this experimental
%   condition. 
%
% Usage:
%   1. opts = BuildFitOpts()
%   2. opts = BuildFitOpts(opts, newopts)
%   3. opts = BuildFitOpts(opts, model, condition, objectives, newopts)
%
% Pass the vectors of models, conditions, and objectives, plus the opts struct
%   to FitObjective (and ObjectiveValue, ObjectiveGradient, etc.) 
%
% Inputs:
%   opts [ options struct scalar ]
%       Existing options struct
%   model [ model struct scalar ]
%       Model to add to fit
%   condition [ condition struct scalar ]
%       Condition based on model to add to fit
%   objectives [ objective function struct scalar or vector ]
%       Objective function(s) based on model+condition to add to fit
%   newopts [ options struct scalar ]
%       New options corresponding to this model+condition+objective to be built
%       into opts. Alternatively, just passing in newopts adds general options to
%       the fit. Allowed options are below.
%
% Outputs:
%   opts [ options struct scalar ]
%       Updated options built from existing opts plus new options from
%       model+condition+objective+newopts
%
% Options:
%  General options for usage 2 newopts
%   .Verbose [ integer scalar {1} ]
%       Larger integers show more debugging output.
%   .ModelsShareParams [ true | {false} ]
%       Whether different models (classes of models) share the
%       same specs for shared/unique params. If true, all param
%       specs for each of k, s, q, h must be unique among all
%       models. Facilitates easy shared params between different
%       models. If false, each model class is treated
%       independently, individual param specs are
%       unambiguous, and each k, s, q, h for each model class
%       has a self contained param spec. Note: false assumes a
%       single model type.
%   .Normalized [ logical scalar {true} ]
%       Indicates if the optimization should be done in log parameters
%       space
%	.UseAdjoint [ logical scalar {true} ]
%       Indicates whether the gradient should be calculated via the
%       adjoint method or the forward method
%   .RelTol [ positive scalar double {1e-6} ]
%       Relative tolerance of the integration
% 	.TolOptim [ positive scalar {1e-5} ]
%       The objective tolerance for optimization. The optimization stops when it is
%       predicted that the objective function cannot be improved more
%       than this in the next iteration.
% 	.Restart [ nonnegative integer scalar {0} ]
%       A scalar integer determining how many times the optimzation
%       should restart once optimization has stopped.
% 	.RestartJump [ handle @(iter,G) returns nonnegative vector nT or
%           scalar | nonnegative vector nT or scalar {0.001} ]
%       This function handle controls the schedule for the noise that
%       will be added to the parameters before each restart. The
%       parameters for the next iteration will be normally distributed
%       in log space with a mean equal to the previous iteration and a
%       standard deviation equal to the value returned by this
%       function. The value returned should usually be a scalar, but it
%       can also be a vector with length equal to the number of active
%       parameters. It can also be numeric, and the noise will be
%       treated as this constant value.
%  	.TerminalObj [ real scalar {-inf} ]
%       Optimization is halted when this objective function value is
%       reached
%   .MaxStepSize [ nonegative scalar {1} ]
%       Scalar fraction indicator of the maximum relative step size
%       that any parameter can take in a single interation
% 	.Algorithm [ string {active-set} ]
%       Option for fmincon. Which optimization algorithm to use
% 	.MaxIter [ postive scalar integer {1000} ]
%       Option for fmincon. Maximum number of iterations allowed before
%       optimization will be terminated.
% 	.MaxFunEvals [ postive scalar integer {5000} ]
%       Option for fmincon. Maximum number of objective function
%       evaluations allowed before optimization will be terminated.
%   .Method [ {'fmincon'} ]
%       The local optimization method. Right now only fmincon is
%       implemented. Others, including stochastic methods may be
%       implemented later.
%   .ComplexStep [ true | {false} ]
%       Whether to use complex step differentiation. Useful for calculations
%       that don't have analytic derivatives or when checking (n+1)-th order
%       terms.
%   .ComputeSensPlusOne [ true | {false} ]
%       Whether to compute the (n+1)-th order sensitivity for an
%       n-th order calculation. Example: NLME requires 1st order
%       sensitivities (gradient) for the objective function calculation and
%       2nd order sensitivities (Hessian) for fitting.
%	.RunParallelExpts [ true | {false} | scalar positive integer ]
%       Whether to run experiments in parallel. If true,
%       initializes a parallel pool using max number of threads
%       available, or if a scalar positive integer is specified,
%       using that many threads up to the max number of threads
%       available. Running experiments in parallel is useful
%       when the fit needs to integrate many
%       covariates/experimental conditions/patients for each
%       step of the optimization.
%   .GlobalOptimization [ logical scalar {false} ]
%       Use global optimization in addition to fmincon
%	.GlobalOpts [ options struct scalar {} ]
%       TODO: API in progress
%
%  Experiment-specific options for usage 3 newopts
% 	.UseParams [ nk x 1 double vector {ones(nk,1)} ]
%       Specifies rate parameters to be fit, with nonzero
%       entries indicating fit, with unique/shared between
%       conditions denoted by different/same values in each row.
%       Default is all rate params fit, shared between
%       conditions.
%	.UseSeeds [ ns x 1 double vector {zeros{ns,1}} ]
%       Specifies seed parameters to be fit, with nonzero
%       entries indicating fit, with unique/shared between
%       conditions denoted by different/same values in each row.
%       Default is no seed params fit.
% 	.UseInputControls [ nq x 1 double vector {zeros(nq,1)} ]
%       Specifies input control parameters to be fit, with nonzero
%       entries indicating fit, with unique/shared between
%       conditions denoted by different/same values in each row.
%       Default is no input control params fit.
%   .UseDoseControls [ nh x 1 double vector {zeros(nh,1)} ]
%       Specifies dose control parameters to be fit, with nonzero
%       entries indicating fit, with unique/shared between
%       conditions denoted by different/same values in each row.
%       Default is no dose control params fit.
%   .StartingParams [ nk x 1 double vector {m.k from parent model} ]
%       Starting rate parameter values for the fit. Overrides
%       values in the parent model.
%   .StartingSeeds [ ns x 1 double vector {condition.s} ]
%       Starting seed parameter values for the fit. Overrides
%       values in the supplied condition.
%   .StartingInputControls [ nq x 1 double vector {condition.q} ]
%       Starting input control parameter values for the fit. Overrides
%       values in the supplied condition.
%   .StartingDoseControls [ nh x 1 double vector {condition.h} ]
%       Starting dose control parameter values for the fit. Overrides
%       values in the supplied condition.
%   .ParamLowerBound [ nk x 1 double vector {zeros(nk,1)} ]
%       Minimum values rate parameters can take in constrained fit.
%   .ParamUpperBound [ nk x 1 double vector {inf(nk,1)} ]
%       Maximum values rate parameters can take in constrained fit.
%   .SeedLowerBound [ ns x 1 double vector {zeros(ns,1)} ]
%   .SeedUpperBound [ ns x 1 double vector {inf(ns,1)} ]
%   .InputControlLowerBound [ nq x 1 double vector {zeros(nq,1)} ]
%   .InputControlUpperBound [ nq x 1 double vector {inf(nq,1)} ]
%   .DoseControlLowerBound [ nh x 1 double vector {zeros(nh,1)} ]
%   .DoseControlUpperBound [ nh x 1 double vector {inf(nh,1)} ]
%   .ObjWeights [ nObj x 1 double vector {ones(nObj,1)} ]
%       Applies a post evaluation weight on each objective function
%       in terms of how much it will contribute to the final objective
%       function value. These values aren't normalized.
%   .AbsTol [ positive double vector {1e-9} ]
%       Absolute tolerance of the integration.
%   .ParamSpec [ [Name] x [Type,LB,UB,Use] table {empty} ]
%       Table specifying all parameter bounds and fit specs. The
%       parameter name is the row and Type, LB, UB, and Use are
%       columns. Alternative way of specifying fit spec by name
%       in a nice way. Note: if .ParamSpec is present, rate params k won't be
%       fit by default - only those specified as fit in .ParamSpec will be fit.
%
% Notes:
% - The model+condition construct may correspond to model/model variants.
% - Models and conditions are identified by name. If multiple conditions are
%   added corresponding to 1 model (such as different covariates/patients/etc.),
%   they'll all refer to the same model if they have the same name.
% - If a condition is added with the same model as an existing condition but with newopts specifying
%   fitting parameters independently, a specification of a dummy model is added.
%   FitObjective makes this model automatically as needed.
% - TODO: unit tests for the different ways of calling BuildFitOpts and
%   ConvertFitOpts
% - Usage 2 provides some sanity checking for input general options. However, you
%   can directly add options with opt.(key) = val to hack something in.
%
% Notes for FitObjective implementation
% - ComponentMap and the other fields in fit have nObjectives rows to be
%   consistent. In FitObjective, it will be useful to use the unique number of
%   conditions as rows for integration (and recombining the objectives). For each
%   condition that requires a dummy model (AddDummyModel will be compacted in the
%   same way as the other fields), a dummy model based on the model ind will be
%   made.

%% Clean up input args
if nargin == 0
    % do nothing
elseif nargin == 2
    newopts = varargin{1};
elseif nargin == 5
    model = varargin{1};
    condition = varargin{2};
    objectives = varargin{3};
    newopts = varargin{4};
else
    error('KroneckerBio:BuildFitOpts:InvalidInputArguments', '2 or 5 input args expected')
end

%% Default args
opts_ = [];

opts_.Verbose = 1;

opts_.ModelsShareParams = false;

opts_.RelTol           = 1e-5;

opts_.Normalized = true;
opts_.UseAdjoint = true;

opts_.Aeq              = [];
opts_.beq              = [];
opts_.TolOptim         = 1e-5;
opts_.Restart          = 0;
opts_.RestartJump      = 0.001;
opts_.TerminalObj      = -inf;

opts_.MaxStepSize      = 1;
opts_.Algorithm        = 'active-set';
opts_.MaxIter          = 1000;
opts_.MaxFunEvals      = 5000;
opts_.Method           = 'fmincon';

opts_.ComplexStep = false;

opts_.ComputeSensPlusOne = false;

opts_.RunParallelExpts = false;

opts_.GlobalOptimization = false;
opts_.GlobalOpts         = [];

% Global fit default options
globalOpts.Algorithm = 'globalsearch';
globalOpts.StartPointsToRun = 'bounds-ineqs';
globalOpts.nStartPoints = 10;
globalOpts.UseParallel = false;
globalOpts.MaxIter = 1000; % for pattern search; fix to better default

opts_.GlobalOpts = globalOpts;

opts_.fit = [];

%% Handle 0 input args setting default general opts
if nargin == 0
    opts = opts_;
    return
end

%% Handle 2 input args setting general opts
if nargin == 2
    if ~isempty(newopts)
        opts = pastestruct(opts, newopts);
    end
    return
end

%% Handle 5 input args adding stuff to fit
% Validate inputs
assert(numel(model) == 1);
assert(numel(condition) == 1);
objectives = vec(objectives);
nNewObjectives = length(objectives);

% Set default newopts for the added model+condition+objective
nk = model.nk;
ns = condition.ns;
nq = condition.nq;
nh = condition.nh;
kNames = vec({model.Parameters.Name});
sNames = vec({model.Seeds.Name});

newopts_ = [];
newopts_.UseParams = ones(nk,1);
newopts_.UseSeeds = zeros(ns,1);
newopts_.UseInputControls = zeros(nq,1);
newopts_.UseDoseControls = zeros(nh,1);
newopts_.StartingParams = model.k;
newopts_.StartingSeeds = condition.s;
newopts_.StartingInputControls = condition.q;
newopts_.StartingDoseControls = condition.h;
newopts_.ParamLowerBound        = zeros(nk,1);
newopts_.ParamUpperBound        = inf(nk,1);
newopts_.SeedLowerBound         = zeros(ns,1);
newopts_.SeedUpperBound         = inf(ns,1);
newopts_.InputControlLowerBound = zeros(nq,1);
newopts_.InputControlUpperBound = inf(nq,1);
newopts_.DoseControlLowerBound  = zeros(nh,1);
newopts_.DoseControlUpperBound  = inf(nh,1);
newopts_.ObjWeights = ones(nNewObjectives,1);
newopts_.AbsTol = 1e-9;

newopts = mergestruct(newopts_, newopts);

%% Handle ParamSpec if specified
if isfield(newopts, 'ParamSpec')
    ps = newopts.ParamSpec;
    np = height(ps);
    names = ps.Properties.RowNames;
    for ip = 1:np
        p = ps(ip,:);
        
        name = names{ip};
        type = p.Type{1};
        lb   = p.LB;
        ub   = p.UB;
        use  = p.Use;
        
        switch type
            case 'k'
                ind = find(ismember(kNames, name));
                newopts.ParamLowerBound(ind) = lb;
                newopts.ParamUpperBound(ind) = ub;
                newopts.UseParams(ind) = use;
            case 's'
                ind = find(ismember(sNames, name));
                newopts.SeedLowerBound(ind) = lb;
                newopts.SeedUpperBound(ind) = ub;
                newopts.UseSeeds(ind) = use;
            case 'q'
                error('Not implemented yet, or input control param names not specified anywhere')
            case 'h'
                error('Not implemented yet, or dose control param names not specified anywhere')
            otherwise
                warning('BuildFitObject:ProcessParamSpec:InvalidParamType', 'Param type %s in opts.ParamSpec not recognized. Skipping.', type)
        end
    end
end


%% Validate inputs
newopts.UseParams = vec(newopts.UseParams);
newopts.UseSeeds = vec(newopts.UseSeeds);
newopts.UseInputControls = vec(newopts.UseInputControls);
newopts.UseDoseControls = vec(newopts.UseDoseControls);
assert(length(newopts.UseParams) == nk);
assert(length(newopts.UseSeeds) == ns);
assert(length(newopts.UseInputControls) == nq);
assert(length(newopts.UseDoseControls) == nh);

ParamSpec_ = cell(1,4);
ParamSpec_{1} = newopts.UseParams;
ParamSpec_{2} = newopts.UseSeeds;
ParamSpec_{3} = newopts.UseInputControls;
ParamSpec_{4} = newopts.UseDoseControls;

LowerBounds_ = cell(1,4);
UpperBounds_ = cell(1,4);
LowerBounds_{1} = fixBounds(newopts.ParamLowerBound, nk);
UpperBounds_{1} = fixBounds(newopts.ParamUpperBound, nk);
LowerBounds_{2} = fixBounds(newopts.SeedLowerBound, ns);
UpperBounds_{2} = fixBounds(newopts.SeedUpperBound, ns);
LowerBounds_{3} = fixBounds(newopts.InputControlLowerBound, nq);
UpperBounds_{3} = fixBounds(newopts.InputControlUpperBound, nq);
LowerBounds_{4} = fixBounds(newopts.DoseControlLowerBound, nh);
UpperBounds_{4} = fixBounds(newopts.DoseControlUpperBound, nh);

ParamStart_ = cell(1,4);
ParamStart_{1} = validateBounds(newopts.StartingParams, LowerBounds_{1}, UpperBounds_{1}, kNames);
ParamStart_{2} = validateBounds(newopts.StartingSeeds, LowerBounds_{2}, UpperBounds_{2}, sNames);
ParamStart_{3} = validateBounds(newopts.StartingInputControls, LowerBounds_{3}, UpperBounds_{3});
ParamStart_{4} = validateBounds(newopts.StartingDoseControls, LowerBounds_{4}, UpperBounds_{4});

newopts.ObjWeights = vec(newopts.ObjWeights);
assert(length(newopts.ObjWeights) == nNewObjectives);

newopts.Continuous = any([objectives.Continuous]);
newopts.Complex = any([objectives.Complex]);
newopts.tGet = [];
for iObj = 1:nNewObjectives
    newopts.tGet = [newopts.tGet; vec(objectives(iObj).DiscreteTimes)];
end
newopts.tGet = unique(newopts.tGet);

if ~isscalar(newopts.AbsTol) || ~isnumeric(newopts.AbsTol)
    warning('KroneckerBio:BuildFitOpts:AbsTolFormNotRecognized', 'BuildFitOpts currently only supports a scalar number of AbsTol. Hope you know what you''re doing.')
end

%% Get private fit struct from opts (holds useful fields)
if isfield(opts, 'fit')
    fit = opts.fit;
else
    fit = [];
end

if isfield(fit, 'ModelNames')
    ModelNames = fit.ModelNames;
    assert(iscellstr(ModelNames))
else
    ModelNames = {};
end
nModels = length(ModelNames);

if isfield(fit, 'ConditionNames')
    ConditionNames = fit.ConditionNames;
    assert(iscellstr(ConditionNames))
else
    ConditionNames = {};
end
nConditions = length(ConditionNames);

if isfield(fit, 'ObjectiveNames')
    ObjectiveNames = fit.ObjectiveNames;
    assert(iscellstr(ObjectiveNames))
else
    ObjectiveNames = {};
end
nObjectives = length(ObjectiveNames);

if isfield(fit, 'UseParams')
    UseParams = fit.UseParams;
    assert(iscell(UseParams) && length(UseParams) == nConditions)
else
    UseParams = {};
end

if isfield(fit, 'UseSeeds')
    UseSeeds = fit.UseSeeds;
    assert(iscell(UseSeeds) && length(UseSeeds) == nConditions)
else
    UseSeeds = {};
end

if isfield(fit, 'UseInputControls')
    UseInputControls = fit.UseInputControls;
    assert(iscell(UseInputControls) && length(UseInputControls) == nConditions)
else
    UseInputControls = {};
end

if isfield(fit, 'UseDoseControls')
    UseDoseControls = fit.UseDoseControls;
    assert(iscell(UseDoseControls) && length(UseDoseControls) == nConditions)
else
    UseDoseControls = {};
end

if isfield(fit, 'ParamSpec')
    ParamSpec = fit.ParamSpec;
    assert(iscell(ParamSpec) && length(ParamSpec) == nConditions)
else
    ParamSpec = {};
end

if isfield(fit, 'ParamStart')
    ParamStart = fit.ParamStart;
    assert(iscell(ParamStart) && length(ParamStart) == nConditions)
else
    ParamStart = {};
end

if isfield(fit, 'LowerBounds')
    LowerBounds = fit.LowerBounds;
    assert(iscell(LowerBounds) && length(LowerBounds) == nConditions)
else
    LowerBounds = {};
end

if isfield(fit, 'UpperBounds')
    UpperBounds = fit.UpperBounds;
    assert(iscell(UpperBounds) && length(UpperBounds) == nConditions)
else
    UpperBounds = {};
end

if isfield(fit, 'Continuous')
    Continuous = fit.Continuous;
    assert(islogical(Continuous) && length(Continuous) == nConditions)
else
    Continuous = false(0,1);
end

if isfield(fit, 'Complex')
    Complex = fit.Complex;
    assert(islogical(Complex) && length(Complex) == nConditions)
else
    Complex = false(0,1);
end

if isfield(fit, 'tGet')
    tGet = fit.tGet;
    assert(iscell(tGet) && length(tGet) == nConditions)
else
    tGet = {};
end

if isfield(fit, 'ObjWeights')
    ObjWeights = fit.ObjWeights;
    assert(iscell(ObjWeights) && length(ObjWeights) == nConditions)
else
    ObjWeights = {};
end

if isfield(fit, 'AbsTol')
    AbsTol = fit.AbsTol;
    assert(iscell(AbsTol) && length(AbsTol) == nConditions)
else
    AbsTol = {};
end

if isfield(fit, 'AddDummyModel')
    AddDummyModel = fit.AddDummyModel;
    assert(isnumeric(AddDummyModel) && length(AddDummyModel) == nConditions)
else
    AddDummyModel = [];
end

%% Get or make new component map (maintains model+condition+objective mapping)
% Each objective will add a row to the component map
if isfield(fit, 'ComponentMap')
    ComponentMap = fit.ComponentMap;
    assert(iscell(ComponentMap) && all(size(ComponentMap) == [nConditions,3]))
else
    ComponentMap = cell(0,3);
end

if isfield(fit, 'ParamMapper')
    ParamMapper = fit.ParamMapper;
    assert(isa(ParamMapper, 'ParamMapper'))
else
    if opts.ModelsShareParams
        ParamMapper = ParamMapperMultiModelType;
    else
        ParamMapper = ParamMapperOneModelType;
    end
end

%% Add objective(s)
newObjectiveNames = vec({objectives.Name});
newObjectiveInds = zeros(nNewObjectives,1);
for i = 1:nNewObjectives
    newObjectiveInds(i) = length(ObjectiveNames) + i;
end

%% Add condition
newConditionName = condition.Name;
newConditionInd = length(ConditionNames) + 1;

%% Add model or see if model already exists
newModelName = model.Name;
existingModelInd = find(ismember(ModelNames, newModelName)); % can be multiple matches if multiple conditions share a model
needDummyModel = true;
if isempty(existingModelInd) % new model
    newModelInd = length(ModelNames) + 1; % assign new model index
    needDummyModel = false;
else % existing model
    % Check if existing model shares UseParams. If existing model shares
    %   UseParams, then the new condition will refer to the existing model. If the
    %   existing model has a different UseParams, make a note for FitObjective
    %   to add a dummy model when used.
    for i = 1:length(existingModelInd)
        existingUseParams = UseParams{existingModelInd};
        if all(newopts.UseParams == existingUseParams)
            newModelInd = existingModelInd(1);
            needDummyModel = false;
        end
    end
    % If it falls through through to here w/o a match, needDummyModel remains true
end

if needDummyModel
    newModelInd = existingModelInd(1);
    newDummyInd = newModelInd;
else
    newDummyInd = 0;
end

%% Assemble fit fields
i = nConditions + 1;
ModelNames{i} = newModelName;
ConditionNames{i} = newConditionName;
UseParams{i} = newopts.UseParams;
UseSeeds{i} = newopts.UseSeeds;
UseInputControls{i} = newopts.UseInputControls;
UseDoseControls{i} = newopts.UseDoseControls;
ParamSpec{i} = ParamSpec_;
ParamStart{i} = ParamStart_;
LowerBounds{i} = LowerBounds_;
UpperBounds{i} = UpperBounds_;
Continuous(i) = newopts.Continuous;
Complex(i) = newopts.Complex;
tGet{i} = newopts.tGet;
AbsTol{i} = newopts.AbsTol;
ObjWeights{i} = newopts.ObjWeights;
AddDummyModel(i) = newDummyInd;
ObjectiveNames = [ObjectiveNames; newObjectiveNames];

%% Assemble component map
ComponentMap = [ComponentMap; {newModelInd, newConditionInd, newObjectiveInds}];

%% Assemble parameter mapper
ParamMapper.AddCondition({newopts.UseParams, newopts.UseSeeds, newopts.UseInputControls, newopts.UseDoseControls});

%% Pack options back into opts
fit.ModelNames = vec(ModelNames);
fit.ConditionNames = vec(ConditionNames);
fit.ObjectiveNames = vec(ObjectiveNames);
fit.UseParams = vec(UseParams);
fit.UseSeeds = vec(UseSeeds);
fit.UseInputControls = vec(UseInputControls);
fit.UseDoseControls = vec(UseDoseControls);
fit.ParamSpec = vec(ParamSpec);
fit.ParamStart = vec(ParamStart);
fit.LowerBounds = vec(LowerBounds);
fit.UpperBounds = vec(UpperBounds);
fit.Continuous = vec(Continuous);
fit.Complex = vec(Complex);
fit.tGet = vec(tGet);
fit.AbsTol = vec(AbsTol);
fit.ObjWeights = vec(ObjWeights);
fit.AddDummyModel = vec(AddDummyModel);
fit.ComponentMap = ComponentMap;
fit.ParamMapper = ParamMapper;
fit.T2Tlocal = @ParamMapper.T2Tlocal;
fit.Tlocal2T = @ParamMapper.Tlocal2T;
opts.fit = fit;
% opts.Split = @split;


% %% Function for returning a new opts struct corresponding to a subset of conditions
%     function opts_i = split(iw, nw)
%         % Split experiment-specific fields in opts (opts.fit) to make a new fitting
%         %   scheme for a subset of the experiments. Useful for splitting experiments
%         %   across workers to be integrated in parallel.
%         %
%         % Inputs:
%         %   iw [ positive scalar integer ]
%         %   nw [ positive scalar integer ]
%         %
%         % Outputs:
%         %   opts_i [ options struct scalar ]
%         componentsMaps = splitComponentMap(opts.fit.ComponentMap, nw);
%         componentMap = componentsMaps{iw};
%         
%         if isempty(componentMap)
%             opts_i = opts;
%             opts_i.fit = [];
%             return
%         end
%         
%         
%     end

end

%% Helper functions
function fixedBounds = fixBounds(bounds, n)
% Convert bounds into standard format.
%
% Inputs:
%   bounds [ double scalar | double vector ]
%       Value for bounds. Scalar is expanded.
%   n [ double scalar ]
%       Target number of bounds
%
% Outputs:
%   fixedBounds [ double vector n ]
%       Fixed (copied, column vector) bounds.
% Throws error if invalid number of bounds supplied.
nBounds = length(bounds);
switch nBounds
    case 1
        fixedBounds = repmat(bounds, n, 1);
    case n
        fixedBounds = vec(bounds);
    otherwise
        error('FitObject:fixBounds: Number of bounds %d not equal to 1 or number params %d', nBounds, n)
end
end

function out = validateBounds(in, lb, ub, names)
% Check bounds on parameters: returning it if it satisfies them, setting equal to
% appropriate bound if it falls outside and throwing warning. Accepts vectors of
% in, lb, and ub. Throws error if input vectors' lengths don't match.
%
% Inputs:
%   in [ nx x 1 double vector ]
%       Vector of parameters
%   lb [ nx x 1 double vector ]
%       Vector of lower bounds
%   ub [ nx x 1 double vector ]
%       Vector of upper bounds
%   names [ nx x 1 cell vector of strings ]
%       Optional names of parameters for useful warning/error messages
%
% Outputs:
%   out [ nx x 1 double vector ]
%       Vector of validated/floored+ceilinged parameters

if nargin < 4
    names = [];
end

% Make sure all inputs have the same dimension
in = vec(in);
lb = vec(lb);
ub = vec(ub);
nin = length(in);
nlb = length(lb);
nub = length(ub);
assert(nin == nlb, 'BuildFitOpts:validateBounds:lbLengthMismatch', 'Fit param vector length %d and lower bounds length %d don''t match', nin, nlb);
assert(nin == nub, 'BuildFitOpts:validateBounds:ubLengthMismatch', 'Fit param vector length %d and upper bounds length %d don''t match', nin, nub);
if ~isempty(names)
    nNames = length(names);
    assert(nin == nNames, 'BuildFitOpts:validateBounds:nameLengthMismatch', 'Fit param vector length %d and names length %d don''t match', nin, nNames);
end

lMask = in < lb;
if any(lMask)
    lInds = find(lMask);
    for i = 1:numel(lInds)
        lInd = lInds(i);
        if isempty(names)
            name = ['param ' num2str(lInd)];
        else
            name = names{lInd};
        end
        warning('BuildFitOpts:validateBounds: fit param %s = %d below lower bound = %d. Setting to lower bound', name, in(lInd), lb(lInd))
    end
end
in(lMask) = lb(lMask);

uMask = in > ub;
if any(uMask)
    uInds = find(uMask);
    for i = 1:numel(uInds)
        uInd = uInds(i);
        if isempty(names)
            name = ['param ' num2str(uInd)];
        else
            name = names{uInd};
        end
        warning('FitObject:validateBounds: fit param %s = %d above upper bound = %d. Setting to upper bound', name, in(uInd), ub(uInd))
    end
end
in(uMask) = ub(uMask);

out = in;

end



