function opts = BuildFitOpts(opts, varargin)
% Build options for advanced fitting schemes. Each set of
%   model+condition+objective corresponds to an experimental condition, it's
%   dependent model, and all the objective functions for this experimental
%   condition. 
%
% Usage:
%   opts = BuildFitOpts(opts, newopts)
%   opts = BuildFitOpts(opts, model, condition, objectives, newopts)
%
% Pass the vectors of models, conditions, and objectives, plus the opts struct
%   to FitObjective (and ObjectiveValue, ObjectiveGradient, etc.) 
%
% Inputs:
%   opts [ options struct scalar ]
%       Existing options struct
%   model [ model struct scalar ]
%       Model 
%   condition [ condition struct scalar ]
%
%   objectives [ objective function struct scalar or vector ]
%
%   newopts [ options struct scalar ]
%       New options corresponding to this model+condition+objective to be built
%       into opts. Alternatively, just passing in newopts adds global options to
%       the fit.
%
% Outputs:
%   opts [ options struct scalar ]
%       Updated options built from existing opts plus new options from
%       model+condition+objective+newopts
%
% Notes:
% - The model+condition construct may correspond to model/model variants.
% - Models and conditions are identified by name. If multiple conditions are
% added corresponding to 1 model (such as different covariates/patients/etc.),
% they'll all refer to the same model if they have the same name.
% - If a condition is added with the same model as an existing condition but with newopts specifying
% fitting parameters independently, a specification of a dummy model is added.
% FitObjective makes this model automatically as needed.
% - TODO: unit tests for the different ways of calling BuildFitOpts and
% ConvertFitOpts
%
% Notes for FitObjective implementation
% - ComponentMap and the other fields in fit have nObjectives rows to be
% consistent. In FitObjective, it will be useful to use the unique number of
% conditions as rows for integration (and recombining the objectives). For each
% condition that requires a dummy model (AddDummyModel will be compacted in the
% same way as the other fields), a dummy model based on the model ind will be
% made.

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

%% Handle 0 input args setting default global opts
if nargin == 0
    opts = [];
    
    opts.Verbose = 1;
    
    opts.RelTol           = 1e-5;
    opts.AbsTol           = [];
    
    opts.Normalized = true;
    opts.UseAdjoint = true;
    
    opts.Aeq              = [];
    opts.beq              = [];
    opts.TolOptim         = 1e-5;
    opts.Restart          = 0;
    opts.RestartJump      = 0.001;
    opts.TerminalObj      = -inf;
    
    opts.MaxStepSize      = 1;
    opts.Algorithm        = 'active-set';
    opts.MaxIter          = 1000;
    opts.MaxFunEvals      = 5000;
    opts.Method           = 'fmincon';
    
    opts.ComputeSensPlusOne = false;
    
    opts.RunParallelExpts = false;
    
    opts.UseIndependentModels = false;
    
    opts.GlobalOptimization = false;
    opts.GlobalOpts         = [];
    
    % Global fit default options
    globalOpts.Algorithm = 'globalsearch';
    globalOpts.StartPointsToRun = 'bounds-ineqs';
    globalOpts.nStartPoints = 10;
    globalOpts.UseParallel = false;
    globalOpts.MaxIter = 1000; % for pattern search; fix to better default
    
    opts.GlobalOpts = globalOpts;
    
    return
end

%% Handle 2 input args setting global opts
if nargin == 2
    if ~isempty(newopts)
        fields = fieldnames(newopts);
        for iField = 1:length(fields)
            field = fields{iField};
            opts.(field) = newopts.(field);
        end
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

newopts_ = [];
newopts_.UseParams = ones(nk,1);
newopts_.UseSeeds = zeros(ns,1);
newopts_.UseInputControls = zeros(nq,1);
newopts_.UseDoseControls = zeros(nh,1);
newopts_.ParamLowerBound        = zeros(nk,1);
newopts_.ParamUpperBound        = inf(nk,1);
newopts_.SeedLowerBound         = zeros(ns,1);
newopts_.SeedUpperBound         = inf(ns,1);
newopts_.InputControlLowerBound = zeros(nq,1);
newopts_.InputControlUpperBound = inf(nq,1);
newopts_.DoseControlLowerBound  = zeros(nh,1);
newopts_.DoseControlUpperBound  = inf(nh,1);
newopts_.ObjWeights = ones(nNewObjectives,1);

newopts = mergestruct(newopts_, newopts);

% Validate inputs
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

newopts.ObjWeights = vec(newopts.ObjWeights);
assert(length(newopts.ObjWeights) == nNewObjectives);

newopts.Continuous = any([objectives.Continuous]);
newopts.Complex = any([objectives.Complex]);
newopts.tGet = [];
for iObj = 1:nNewObjectives
    newopts.tGet = [newopts.tGet; vec(objectives(iObj).DiscreteTimes)];
end
newopts.tGet = unique(newopts.tGet);

if isfield(newopts, 'AbsTol')
    newopts.AbsTol = fixAbsTol(newopts.AbsTol, 2, newopts.continuous, model.nx, 1, opts.UseAdjoint, logical(newopts.UseParams), logical(newopts.UseSeeds), logical(newopts.UseInputControls), logical(newopts.UseDoseControls));
else
    newopts.AbsTol = 1e-9;
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

if isfield(fit, 'AbsTol')
    AbsTol = fit.AbsTol;
    assert(iscell(AbsTol) && length(AbsTol) == nConditions)
else
    AbsTol = {};
end

if isfield(fit, 'ObjWeights')
    ObjWeights = fit.ObjWeights;
    assert(iscell(ObjWeights) && length(ObjWeights) == nConditions)
else
    ObjWeights = {};
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
    assert(isnumeric(ComponentMap) && all(size(ComponentMap) == [nObjectives,3]))
else
    ComponentMap = zeros(0,3);
end

if isfield(fit, 'ParamMapper')
    ParamMapper = fit.ParamMapper;
    assert(isa(ParamMapper, 'ParamMapper'))
else
    if opts.UseIndependentModels
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
        if newopts.UseParams == existingUseParams
            newModelInd = existingModelInd;
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
for i = 1:nNewObjectives
    ComponentMap = [ComponentMap; newModelInd, newConditionInd, newObjectiveInds(i)];
end

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





