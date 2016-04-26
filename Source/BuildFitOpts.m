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
%
% Notes for FitObjective implementation
% - ComponentMap and the other fields in fit have nObjectives rows to be
% consistent. In FitObjective, it will be useful to use the unique number of
% conditions as rows for integration (and recombining the objectives). For each
% condition that requires a dummy model (AddDummyModel will be compacted in the
% same way as the other fields), a dummy model based on the model ind will be
% made.

%% Clean up input args
if nargin == 2
    newopts = varargin{1};
elseif nargin == 5
    model = varargin{1};
    condition = varargin{2};
    objectives = varargin{3};
    newopts = varargin{4};
else
    error('KroneckerBio:BuildFitOpts:InvalidInputArguments', '2 or 5 input args expected')
end

%% Handle 2 input args setting global opts
if nargin == 2
    fields = fieldnames(newopts);
    for iField = 1:length(fields)
        field = fields{iField};
        opts.(field) = newopts.(field);
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
newopts_.UseInputControls = zeros(nq, 1);
newopts_.UseDoseControls = zeros(nh, 1);

newopts = mergestruct(newopts_, newopts);

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
    assert(iscell(UseParams) && length(UseParams) == nObjectives)
else
    UseParams = {};
end

if isfield(fit, 'UseSeeds')
    UseSeeds = fit.UseSeeds;
    assert(iscell(UseSeeds) && length(UseSeeds) == nObjectives)
else
    UseSeeds = {};
end

if isfield(fit, 'UseInputControls')
    UseInputControls = fit.UseInputControls;
    assert(iscell(UseInputControls) && length(UseInputControls) == nObjectives)
else
    UseInputControls = {};
end

if isfield(fit, 'UseDoseControls')
    UseDoseControls = fit.UseDoseControls;
    assert(iscell(UseDoseControls) && length(UseDoseControls) == nObjectives)
else
    UseDoseControls = {};
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
for i = 1:nNewObjectives
    ModelNames = [ModelNames; newModelName];
    ConditionNames = [ConditionNames; newConditionName];
    UseParams = [UseParams; newopts.UseParams];
    UseSeeds = [UseSeeds; newopts.UseSeeds];
    UseInputControls = [UseInputControls; newopts.UseInputControls];
    UseDoseControls = [UseDoseControls; newopts.UseDoseControls];
    AddDummyModel = [AddDummyModel; newDummyInd];
end
ObjectiveNames = [ObjectiveNames; newObjectiveNames];


%% Assemble component map
for i = 1:nNewObjectives
    ComponentMap = [ComponentMap; newModelInd, newConditionInd, newObjectiveInds(i)];
end

%% Pack options back into opts
fit.ModelNames = ModelNames;
fit.ConditionNames = ConditionNames;
fit.ObjectiveNames = ObjectiveNames;
fit.UseParams = UseParams;
fit.UseSeeds = UseSeeds;
fit.UseInputControls = UseInputControls;
fit.UseDoseControls = UseDoseControls;
fit.AddDummyModel = AddDummyModel;
fit.ComponentMap = ComponentMap;
opts.fit = fit;
end

%% Helper functions
% Build T2Tlocal and Tlocal2T mappers






