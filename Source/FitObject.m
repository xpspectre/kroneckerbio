classdef FitObject < handle
    % Specifies system for FitObjective
    % Models, Conditions, and Objectives (objective functions) are stored in
    % struct arrays. They can be directly indexed by calling code for
    % convenience.
    % The 1st of each component is the default; other components that are added
    % refer to this one if not otherwise specified.
    % Makes sure entire system is consistent with each component added.
    % Contains all the default options required for fitting (previously in FitObjective)
    
    % TODO/Notes
    % Use property listeners - as soon as something is changed/added, update
    % meta properties like paramsMap
    % Unify param types k, s, q, h. One strategy: make all calls to them
    %   analogous; then combine and eliminate redundancy
    % Having models, conditions, objectives be objects would allow them to be
    %   handle objects/by reference, potentially saving a lot of space with big
    %   models or model variants
    % Some potentially relevant design patterns, esp if models, conditions,
    % objectives are objects:
    %   composite: implement traversible hierarchy of models < conditions < objectives
    %   visitor: implementing additional actions for hierarchy components - may
    %       be excessive
    %   decorator: adding fitting-specific functionality to components w/o
    %       polluting them when they're simply used for simulating
    
   properties
       Name
       modelNames = {}
       conditionNames = {}
       objectiveNames = {}
       nModels = 0
       nConditions = 0;
       nObjectives = 0
       options
       componentMap = zeros(0,3) % nObjectives x 3 double matrix of [modelIdx, conditionIdx, objectiveIdx]; Note: when parts are turned into objects, get rid of this
       paramsMap = cell(0,4) % nObjectives x 4 cell matrix of [param, seed, input, dose] to T index mapping; 0 indicates not fit; multiple refs to a single T index are summed
       paramsShared = cell(0,4) % nObjectives x 4 cell matrix of [param, seed, input, dose] shared params index spec
       modelClassMap = [] % nObjectives x 1 double matrix of model class index each objective belongs to
       nT
   end
   
   properties (SetObservable)
       Models
       Conditions
       Objectives
   end
   
   methods
       function this = FitObject(name, opts)
           % Constructor; Instantiate FitObject.
           % Inputs:
           %    opts [ struct ]
           %        Options struct allowing the following fields:
           %        .Verbose [ intger scalar {tbd} ]
           %            Larger integers show more debugging output.
           %        .ModelsShareParams [ boolean scalar {false} ]
           %            Whether different models (classes of models) share the
           %            same specs for shared/unique params. If true, all param
           %            specs for each of k, s, q, h must be unique among all
           %            models. Facilitates easy shared params between different
           %            models. If false, each model class is treated
           %            independently, individual param specs are
           %            unambiguous, and each k, s, q, h for each model class
           %            has a self contained param spec.
           %            TODO: implement this
           
           % Clean up inputs
           if nargin < 2
               opts = [];
               if nargin < 1
                   name = [];
               end
           end
           
           % Fill in defaults
           if isempty(name)
               name = 'DefaultFitObject';
           end
           
           % Default options
           opts_.Verbose = 2;
           opts_.ModelsShareParams = false;
           
           opts_.Normalized       = true;
           opts_.UseAdjoint       = false; % set to true when Adjoit calculation is updated
           
           opts_.TolOptim         = 1e-5;
           opts_.Restart          = 0;
           opts_.RestartJump      = 0.001;
           opts_.TerminalObj      = -inf;
           
           opts_.MaxStepSize      = 1;
           opts_.Algorithm        = 'active-set';
           opts_.MaxIter          = 1000;
           opts_.MaxFunEvals      = 5000;
           
           opts_.GlobalOptimization = false;
           opts_.GlobalOpts         = [];
           
           opts_.Aeq = [];
           opts_.beq = [];
           
           optsGlobal_.Algorithm = 'globalsearch';
           optsGlobal_.StartPointsToRun = 'bounds-ineqs';
           optsGlobal_.nStartPoints = 10;
           optsGlobal_.UseParallel = false;
           optsGlobal_.MaxIter = 1000; % for pattern search; fix to better default
           
           opts = mergestruct(opts_, opts);
           opts.GlobalOpts = mergestruct(optsGlobal_, opts.GlobalOpts);
           
           % Input validation and checking
           % TODO
           
           % Attach listeners
           this.addlistener('Models', 'PostSet', @(src,event)this.updateModels(src,event));
           this.addlistener('Conditions', 'PostSet', @(src,event)this.updateConditions(src,event));
           this.addlistener('Objectives', 'PostSet', @(src,event)this.updateObjectives(src,event));
           
           % Assign fields
           this.Name = name;
           this.options = opts;
           
       end
       
       function addOptions(this, opts)
           % Add options to fit object
           % TODO: validate inputs
           this.options = mergestruct(this.options, opts);
       end
       
       function addModel(this, model, opts)
           if nargin < 3
               opts = [];
           end
           
           % Default options
           opts_.BaseModel = model.Name; % default is this model is a new type of model
           opts = mergestruct(opts_, opts);
           
           % Validate and annotate models with new information
           model.BaseModel = opts.BaseModel;
           
           % Add model to fit object
           this.Models = [this.Models; model];
       end
       
       function addCondition(this, condition, opts)
           % Add experimental condition to fit object.
           % Automatically associates with model with corresponding name (Note: when turned into objects, enforce this with typing)
           
           if nargin < 3
               opts = [];
           end
           
           % Check that associated model is already in fit object
           parentModelMask = ismember(this.modelNames, condition.ParentModelName);
           if ~parentModelMask
               error('FitObject:addCondition: Condition parent model %s not present in fit object.', condition.ParentModelName)
           end
           
           % Add associated model index to condition for convenience
           parentModelIdx = find(parentModelMask);
           condition.ParentModelIdx = parentModelIdx;
           
           % Correct count of each param type
           nk = this.Models(condition.ParentModelIdx).nk;
           ns = condition.ns;
           nq = condition.nq;
           nh = condition.nh;
           
           % Existing params
           kStart = this.Models(condition.ParentModelIdx).k;
           sStart = condition.s;
           qStart = condition.q;
           hStart = condition.h;
           
           % Default options
           opts_ = [];
           opts_.UseParams = ones(nk, 1); % fit all rate params, same in model class
           opts_.UseSeeds = ones(ns, 1); % fit all seeds, same in model class
           opts_.UseInputControls = zeros(nq, 1); % don't fit input control params
           opts_.UseDoseControls = zeros(nh, 1); % don't fit dose control params
           
           opts_.StartingParams = [];
           opts_.StartingSeeds = [];
           opts_.StartingInputControls = [];
           opts_.StartingDoseControls = [];
           
           opts_.ParamLowerBound = 0;
           opts_.ParamUpperBound = Inf;
           opts_.SeedLowerBound = 0;
           opts_.SeedUpperBound = Inf;
           opts_.InputControlLowerBound = 0;
           opts_.InputControlUpperBound = Inf;
           opts_.DoseControlLowerBound = 0;
           opts_.DoseControlUpperBound = Inf;
           
           opts_.RelTol = 1e-6;
           opts_.AbsTol = 1e-9;
           
           opts = mergestruct(opts_, opts);
           
           % Validate bound and add to condition
           bounds = cell(2,4);
           bounds{1,1} = fixBounds(opts.ParamLowerBound, nk);
           bounds{2,1} = fixBounds(opts.ParamUpperBound, nk);
           bounds{1,2} = fixBounds(opts.SeedLowerBound, ns);
           bounds{2,2} = fixBounds(opts.SeedUpperBound, ns);
           bounds{1,3} = fixBounds(opts.InputControlLowerBound, nq);
           bounds{2,3} = fixBounds(opts.InputControlUpperBound, nq);
           bounds{1,4} = fixBounds(opts.DoseControlLowerBound, nh);
           bounds{2,4} = fixBounds(opts.DoseControlUpperBound, nh);
           condition.Bounds = bounds;
           
           % Sanity check: make sure useX are valid for model and condition
           assert(nk == length(opts.UseParams), 'FitObject:addCondition: Number of UseParams doesn''t match number of model rate params')
           assert(ns == length(opts.UseSeeds), 'FitObject:addCondition: Number of UseSeeds doesn''t match number of model seed params')
           assert(nq == length(opts.UseInputControls), 'FitObject:addCondition: Number of UseInputControls doesn''t match number of model input control params')
           assert(nh == length(opts.UseDoseControls), 'FitObject:addCondition: Number of UseDoseControls doesn''t match number of model dose control params')
           
           % Validate params specs and add to condition
           condition.ParamsSpec = {fixParamsSpec(opts.UseParams, nk), ...
               fixParamsSpec(opts.UseSeeds, ns), ...
               fixParamsSpec(opts.UseInputControls, nq), ...
               fixParamsSpec(opts.UseDoseControls, nh)};
           
           % Override starting params that are fit with user specified values
           % TODO: validate that starting condition falls between bounds; throw
           % warning if it's changed to be within the bounds
           kUse = logical(opts.UseParams);
           sUse = logical(opts.UseSeeds);
           qUse = logical(opts.UseInputControls);
           hUse = logical(opts.UseDoseControls);
           
           if ~isempty(opts.StartingParams) && any(kUse)
               kStart(kUse) = opts.StartingParams(kUse);
           end
           if ~isempty(opts.StartingSeeds) && any(sUse)
               sStart(sUse) = opts.StartingSeeds(sUse);
           end
           if ~isempty(opts.UseInputControls) && any(qUse)
               qStart(qUse) = opts.UseInputControls(qUse);
           end
           if ~isempty(opts.StartingDoseControls) && any(hUse)
               hStart(hUse) = opts.StartingDoseControls(hUse);
           end
           
           if ~isempty(opts.StartingParams) && any(kUse)
               modelName = this.Models(parentModelIdx).Name;
               this.Models(parentModelIdx) = mergestruct(this.Models(parentModelIdx), this.Models(parentModelIdx).Update(kStart));
               this.Models(parentModelIdx).Name = modelName; % hack needed because externally modified fields are lost on model Update
           end
           if (~isempty(opts.StartingSeeds) && any(sUse)) || (~isempty(opts.UseInputControls) && any(qUse)) || (~isempty(opts.StartingDoseControls) && any(hUse))
               condition = mergestruct(condition, condition.Update(sStart, qStart, hStart));
           end
           
           % Add tolerances
           condition.RelTol = opts_.RelTol;
           condition.AbsTol = opts_.AbsTol;
           
           % Add condition to fit object
           this.Conditions = [this.Conditions; condition];
           
       end
       
       function addObjective(this, objective, parentConditions, opts)
           % Add dataset/objective function to fit object.
           % Inputs:
           %    objective [ objective struct ]
           %    parentConditions [ string | cell array of strings {1st condition} ]
           %        names of conditions to associate with
           %    opts [ struct ]
           %        .Weight [ double {1} ]
           %            Relative weight of this objective function when multiple
           %            are summed together in fit object. Simply acts as a
           %            multiplier (i.e., no normalization is needed or added)
           
           if nargin < 4
               opts = [];
               if nargin < 3
                   parentConditions = [];
               end
           end
           
           opts_.Weight = 1;
           opts = mergestruct(opts_, opts);
           
           if isempty(parentConditions)
               parentConditions = {this.Conditions(1).Name};
           elseif ischar(parentConditions)
               parentConditions = {parentConditions};
           end
           
           % Sanity check: make sure parent conditions exist - warn if not
           % found and exit if none are found
           missingConditions = setdiff(parentConditions, this.conditionNames);
           if ~isempty(missingConditions)
               warning('FitObject:addObjective: Conditions %s not found in fit object', cellstr2str(missingConditions))
               
               nMissingConditions = length(missingConditions);
               if nMissingConditions == length(parentConditions)
                   error('FitObject:addObjective: No valid conditions specified by objective function')
               end
           end
           
           % Copy objective function for each condition that it applies to
           conditionsIdxs = find(ismember(this.conditionNames, parentConditions));
           nAddedConditions = length(conditionsIdxs);
           for i = 1:nAddedConditions
               objective_ = objective;
               
               % Annotate objective with weight
               objective_.Weight = opts.Weight;
               
               % Annotate objective with parent condition it's applied to
               parentConditionName = this.conditionNames{conditionsIdxs(i)};
               objective_.ParentConditionName = parentConditionName;
               
               % Sanity check: make sure output indices in objective map to actual
               % outputs in model (referenced thru condition)
               parentConditionIdx = ismember(this.conditionNames, parentConditionName);
               parentModelIdx = ismember(this.modelNames, this.Conditions(parentConditionIdx).ParentModelName);
               parentModelName = this.Models(parentModelIdx).Name;
               outputNames = {this.Models(parentModelIdx).Outputs.Name};
               try
                   fitOutputs = outputNames(objective_.Outputs);
               catch
                   error('FitObject:addObjective: Invalid outputs entered')
               end
               
               % Add objective to fit object
               this.Objectives = [this.Objectives; objective_];
               
               % If verbose, print brief description of outputs being fit.
               if this.options.Verbose >= 2
                   fprintf('Added objective %s that fits outputs %s in condition %s in model %s\n', objective_.Name, cellstr2str(fitOutputs), parentConditionName, parentModelName)
               end
           end
           
       end
       
       function addFitConditionData(this, objective, condition, opts)
           % Higher level/convenience function that adds a "participant" to fit. Abstracts
           % away some of the built in structure that kroneckerbio uses.
           % Inputs:
           %    objective: dataset/objective function of participant
           %    condition: experimental condition, with seeds, inputs, and doses
           % 	opts: options, including indexes or names of each parameter
           %        types' shared/individual grouping
           %        .UseParams [ double vector nk ]
           %            Different numbers indicate independent params. 0
           %            indicates not fit. 
           % Side Effects:
           %    Adds a dummy model if params are allowed to vary between
           %    conditions.
           %
           % Selects 1st condition (and it's corresponding model) if not
           % specified.
           % Only UseParams (if fitting params) needs to be specified in order
           % to see if a dummy model needs to be added.
           % TODO: consider changing this to allow param specs to apply
           % globally. This change may make implementation easier, at the
           % expense of making some simple cases harder for the user. This
           % change would also obviate the need to linear equality constraints.
           if nargin < 4
               opts = [];
               if nargin < 3
                   condition = [];
               end
           end
           
           if isempty(condition)
               condition = this.Conditions(1);
           end
           
           % If params are specified to be fit, verify and possibly add dummy model
           if isfield(opts, 'UseParams')
               UseParams = vec(opts.UseParams);
               assert(this.Models(ismember(this.modelNames, condition.ParentModelName)).nk == length(UseParams), 'FitObject:addFitConditionData: Number of UseParams doesn''t match number of model rate params')
               
               % Add dummy model if k is to be fit independently
               modelIdx = this.addModelOnUniqueParam(UseParams, condition.ParentModelName);
               condition.ParentModelName = this.modelNames{modelIdx};
               
               opts.UseParams = UseParams;
           end
           
           % Add components to fit
           this.addCondition(condition, opts);
           this.addObjective(objective, condition.Name, opts);
       end
       
       function modelIdx = addModelOnUniqueParam(this, useParams, parentModelName)
           % Add a new "dummy" model if useParams specifies unique rate parameters
           % from any existing for the specified parent modelName.
           % Returns index of model that useParams will refer to, either for an
           % existing model if they match or a new one if it needed to be
           % generated.
           % Inputs:
           %    useParams [ nk x 1 double vector ]
           %        Shared/unique param spec
           %    parentModelName [ string ]
           %        String name of parent model class this param set refers to
           % Outputs:
           %    modelIdx [ scalar integer ]
           %        Index of the model useParams/modelName should refer to
           % Side Effects:
           %    Adds a new model to the fit object if needed (trying to fit any unique k)
               
           % Get base model class (index) of specified modelName
           %    Parent model's class will be base model class since base/derived
           %    model hierarchy is limited to depth = 2
           parentModelClass = find(ismember(this.modelNames, parentModelName));
           
           % Get indices of members of this model class
           parentModelClassMembers = find(this.modelClassMap == parentModelClass);
           parentModelClassMembers(parentModelClassMembers==0) = [];
           
           % Get paramsShared of all models in the base model class
           paramsSharedModel = this.paramsShared(this.modelClassMap == parentModelClass);
           
           % First look for exactly matching params specs
           for i = 1:length(paramsSharedModel)
               if all(paramsSharedModel{i} == useParams)
                   modelIdx = parentModelClassMembers(i);
                   return
               end
           end
           
           % Next check if base model has an associated experiment
           % If not, attach to that
           if isempty(ismember(this.componentMap(:,1), parentModelClass))
               modelIdx = parentModelClass;
               return
           end
           
           % No matches - make new model
           model_ = this.getModel(parentModelName);
           model_.Name = [model_.Name '_der' num2str(this.nModels)];
           opts = [];
           opts.BaseModel = parentModelName;
           this.addModel(model_, opts)
           modelIdx = this.nModels; % will refer to index of last model
       end
       
       function model = getModel(this, name)
           % Return model by name.
           modelIdx = ismember(this.modelNames, name);
           if ~any(modelIdx)
               error('FitObject:getModel: Model %s not found', name)
           end
           model = this.Models(modelIdx);
       end
       
       function condition = getCondition(this, name)
           % Return condition by name.
           conditionIdx = ismember(this.conditionNames, name);
           if ~any(conditionIdx)
               error('FitObject:getCondition: Condition %s not found', name)
           end
           condition = this.Conditions(conditionIdx);
       end
       
       
       %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Called inside optimizer's loop
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function [G, D] = computeObjective(this)
           
           % Bridge fit object with computeObj and computeObjGrad forms
           %    Make vector of objectives for each condition
           
           if nargout == 1
               
               G = 0;
               
               for i = 1:this.nConditions
                   modelIdx = this.Conditions(i).ParentModelIdx;
                   objectives = this.Objectives(this.componentMap(:,2) == i)';
                   
                   G = G + computeObj(this.Models(modelIdx), this.Conditions(i), objectives, this.options);
               end
               
           end
           if nargout == 2
               
               G = 0;
               D = zeros(this.nT,1);
               
               for i = 1:this.nConditions
                   modelIdx = this.Conditions(i).ParentModelIdx;
                   objectives = this.Objectives(this.componentMap(:,2) == i)';
                   
                   [Gi, Di] = computeObjGrad(this.Models(modelIdx), this.Conditions(i), objectives, this.options);
                   
                   G = G + Gi;
                   D = mapTlocal2T(Di, this.paramsMap(i,:), this.nT, D, 'add');
               end
               
               % Normalize gradient
               if this.options.Normalized
                   D = D .* this.collectParams;
               end
               
           end
       end
       
       function T = collectParams(this)
           % Collect all params from fit object into master T
           T = zeros(this.nT, 1);
           
           % Go thru all conditions
           for i = 1:this.nConditions
               parentModelIdx = this.Conditions(i).ParentModelIdx;
               Tlocal = {this.Models(parentModelIdx).k, ...
                   this.Conditions(i).s, ...
                   this.Conditions(i).q, ...
                   this.Conditions(i).h};
               T = mapTlocal2T(Tlocal, this.paramsMap(i,:), this.nT, T, 'override');
           end
       end
       
       function updateParams(this, T)
           % Update all params in fit object with master T
           
           % Go thru all conditions
           for i = 1:this.nConditions
               parentModelIdx = this.Conditions(i).ParentModelIdx;
               
               % Get current param values so new values can be pasted over them
               Tlocalold = {this.Models(parentModelIdx).k, this.Conditions(i).s, this.Conditions(i).q, this.Conditions(i).h};
               Tlocal = mapT2Tlocal(T, this.paramsMap(i,:), Tlocalold);
               [k,s,q,h] = deal(Tlocal{:});
               
               
               % mergestruct is exactly what's needed here, but I think this is slow
               % This is a prime reason why we should convert components to objects
               modelName = this.Models(parentModelIdx).Name;
               this.Models(parentModelIdx) = mergestruct(this.Models(parentModelIdx), this.Models(parentModelIdx).Update(k));
               this.Models(parentModelIdx).Name = modelName; % hack needed because externally modified fields are lost on model Update
               this.Conditions(i) = mergestruct(this.Conditions(i), this.Conditions(i).Update(s, q, h));
           end
       end
       
       function [LowerBound, UpperBound] = collectBounds(this)
           % Collect all bounds from fit object into master list
           % Used by normalize and for fmincon options
           LowerBound = zeros(this.nT, 1);
           UpperBound = zeros(this.nT, 1);
           
           % Go thru all conditions
           for i = 1:this.nConditions
               LowerBound = mapTlocal2T(this.Conditions(i).Bounds(1,:), this.paramsMap(i,:), this.nT, LowerBound, 'override');
               UpperBound = mapTlocal2T(this.Conditions(i).Bounds(2,:), this.paramsMap(i,:), this.nT, UpperBound, 'override');
           end
       end
       
       %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Event listeners
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function updateModels(this, ~, ~)
           nModels_ = length(this.Models);
           this.nModels = nModels_;
           
           modelNames_ = cell(nModels_,1);
           for i = 1:nModels_
               modelNames_{i} = this.Models(i).Name;
           end
           this.modelNames = modelNames_;
       end
       
       function updateConditions(this, ~, ~)
           nConditions_ = length(this.Conditions);
           this.nConditions = nConditions_;
           
           conditionNames_ = cell(nConditions_,1);
           for i = 1:nConditions_
               conditionNames_{i} = this.Conditions(i).Name;
           end
           this.conditionNames = conditionNames_;
       end
       
       function updateObjectives(this, ~, ~)
           % Update objectives and all mapping fields for traversing components
           nObjectives_ = length(this.Objectives);
           this.nObjectives = nObjectives_;
           
           objectiveNames_ = cell(nObjectives_,1);
           for i = 1:nObjectives_
               objectiveNames_{i} = this.Objectives(i).Name;
           end
           this.objectiveNames = objectiveNames_;
           
           % Update component map
           this.updateComponentMap;
           
           % Update mapping of objectives to their base model classes
           this.updateModelClassMap;
           
           % Update central params spec locatoin
           this.updateParamsShared;
           
            % Update params mapping and total params vector length
           this.updateParamsMap;
       end
       
       function updateComponentMap(this)
           % Update models < conditions < objectives componentMap which represents the
           % tree as a matrix
           componentMap_ = zeros(this.nObjectives, 3);
           for i = 1:this.nObjectives
               conditionName = this.Objectives(i).ParentConditionName;
               conditionIdx = find(ismember(this.conditionNames, conditionName));
               modelName = this.Conditions(conditionIdx).ParentModelName;
               modelIdx = find(ismember(this.modelNames, modelName));
               componentMap_(i,:) = [modelIdx, conditionIdx, i];
           end
           this.componentMap = componentMap_;
       end
       
       function updateModelClassMap(this)
           % Update mapping of objectives to their base model classes
           modelClassMap_ = zeros(this.nObjectives, 1);
           baseModels = cell(0,1);
           for i = 1:this.nObjectives
               baseModel = this.Models(this.componentMap(i,1)).BaseModel;
               baseModelMask = ismember(baseModels, baseModel);
               if any(baseModelMask)
                   modelClassMap_(i) = find(baseModelMask, 1);
               else
                   baseModels = [baseModels; baseModel];
                   modelClassMap_(i) = max(modelClassMap_) + 1;
               end
           end
           this.modelClassMap = modelClassMap_;
       end
       
       function updateParamsShared(this)
           paramsShared_ = cell(this.nConditions, 4);
           for i = 1:this.nConditions
               paramsShared_(i,:) = this.Conditions(i).ParamsSpec;
           end
           this.paramsShared = paramsShared_;
       end
       
       function updateParamsMap(this)
           [this.paramsMap, this.nT] = getParamsMap(this.paramsShared, this.modelClassMap);
       end
       
   end
   
   methods (Static)
       function fit = buildFitObject(models, conditions, objectives, opts)
           % Factory function to convert legacy format to fit object
           % Inputs:
           %    models [ struct vector nModels ]
           %        Models to fit. Note: ignores all models except 1st one
           %        (since multi-model fitting wasn't previously implemented).
           %    conditions [ struct vector nConditions ]
           %    objectives [ struct matrix nConditions x nObjectives ]
           %    opts [ struct ]
           
           if nargin < 4
               opts = [];
           end
           
           fit = FitObject('Converted_Fit');
           
           % Add models
           nModels_ = length(models);
           if nModels_ > 1
               warning('FitObject:buildFitObject: %d models added - ignoring all but the 1st one', nModels_)
           end
           fit.addModel(models(1));
           
           % Add conditions
           % Assume all conditions belong to 1st model 
           %    (which is necessary because they must all refer to it when being built)
           nConditions_ = length(conditions);
           for i = 1:nConditions_
               fit.addCondition(conditions(i))
           end
           
           % Add objectives
           % Each row of objective functions belongs to the corresponding
           % condition
           nObjectives_ = size(objectives,2);
           for i = 1:nConditions_
               conditionName = conditions(i).Name;
               for j = 1:nObjectives_
                   obj = objectives(i,j);
                   if ~strcmpi(obj.Type, 'Objective.Data') % only add non-objectiveZeros, assuming this generic name isn't overridden
                       fit.addObjective(obj, conditionName);
                   end
               end
           end
           
           % Add options
           fit.addOptions(opts);
       end
   end
   
   
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [paramsMap, nT] = getParamsMap(paramsShared, modelClassMap)
% Update mapping from individual objectives to overall T
Toffset = 0;

kShared = paramsShared(:,1);
kMap = multiModelParamsShared2Map(kShared, modelClassMap);
[kMap, Toffset] = addOffsetToNonzeroCellArray(kMap, Toffset);

sShared = paramsShared(:,2);
sMap = multiModelParamsShared2Map(sShared, modelClassMap);
[sMap, Toffset] = addOffsetToNonzeroCellArray(sMap, Toffset);

qShared = paramsShared(:,3);
qMap = multiModelParamsShared2Map(qShared, modelClassMap);
[qMap, Toffset] = addOffsetToNonzeroCellArray(qMap, Toffset);

hShared = paramsShared(:,4);
hMap = multiModelParamsShared2Map(hShared, modelClassMap);
[hMap, nT] = addOffsetToNonzeroCellArray(hMap, Toffset);

% Updated paramsMap according to (potentially rearranged) specified parameters
paramsMap = [kMap, sMap, qMap, hMap];
end

function [cellOut, offsetOut] = addOffsetToNonzeroCellArray(cellIn, offset)
% Add offset to all nonzero values in cell array
% Also return new offset value
% Assumes cell vertical cell vector
nCells = length(cellIn);
cellOut = cell(nCells,1);
offsetAccumulate = 0;
for i = 1:nCells
    if ~isempty(cellIn{i}) && max(cellIn{i}) > offsetAccumulate
        offsetAccumulate = max(cellIn{i});
    end
    cellOut{i} = (cellIn{i} + offset) .* logical(cellIn{i}); % selecting with logical restores 0's filled in by offset
end
offsetOut = offset + offsetAccumulate;
end

function paramsMap = multiModelParamsShared2Map(paramsShared, modelClassMap)
% Inputs:
%    paramsShared [ cell vector nObjective x 1 ]
%    modelClassMap [ double vector nObjective x 1 ]
%       Indices refer to which groups of models use the same base model
% Outputs:
%    paramsMap [ cell vector nObjective x 1 ]
%
% Note: logically, rate params are treated differently in
% different models while seeds, input control params, and dose
% control params may be more conveniently treated the same in the
% same conditions, regardless of dummy underlying model for
% changed k

% Treat different base models separately
modelClasses = unique(modelClassMap, 'stable');
Toffset = 0;
paramsMap = cell(size(paramsShared));
for modelClass = modelClasses
    
    paramsSharedModel = paramsShared(modelClassMap == modelClass);
    nCons = length(paramsSharedModel);
    nParams = length(paramsSharedModel{1});
    
    paramsSharedMatrix = zeros(nParams, nCons);
    for i = 1:nCons
        paramsSharedMatrix(:,i) = paramsSharedModel{i};
    end
    
    paramsMapMatrix = paramsShared2Map(paramsSharedMatrix) + Toffset;
    paramsMapMatrix(~logical(paramsSharedMatrix)) = 0; % restore 0's to spots overridden by Toffset
    
    paramsMapModel = cell(nCons, 1);
    for i = 1:nCons
        paramsMapModel(i) = {paramsMapMatrix(:,i)};
    end
    
    paramsMap(modelClassMap == modelClass) = paramsMapModel;
    Toffset = max(paramsMapMatrix);
    
end
end

function paramsMap = paramsShared2Map(paramsShared)
% Input:
%    paramsShared: double matrix nk x ncon with rows of
%        shared/unique k
% Output:
%    paramsMap: double matrix nk x ncon with rows of
%        positions in T, relative to 1st fit k in paramsShared
[nk, nCon] = size(paramsShared);
paramsMap = zeros(nk, nCon);
Toffset = 0;
for i = 1:nk
    kUse = paramsShared(i,:);
    kNonZero = logical(kUse);
    kUse = kUse(kNonZero);
    [c,~,ic] = unique(kUse, 'stable');
    
    Tidx = ic' + Toffset;
    
    paramsMap(i,kNonZero) = Tidx;
    Toffset = Toffset + nnz(c);
end
end

function fixedBounds = fixBounds(bounds, n)
% Convert bounds into standard format.
% Inputs:
%   bounds [ double scalar | double vector ]
%       Value for bounds. Scalar is expanded.
%   n [ double scalar ]
%       Target number of bounds
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

function fixedParamsSpec = fixParamsSpec(ParamsSpec, n)
nParams = length(ParamsSpec);
switch nParams
    case n
        fixedParamsSpec = vec(ParamsSpec);
    otherwise
        error('FitObject:fixParamsSpec: Number of params specs %d not equal to number params %d', nParams, n)
end
end

function Tlocal = mapT2Tlocal(T, Tmap, Tlocalold)
% Map overall T vector to local T's in conditions.
% Inputs:
%    T [ nT x 1 double vector ]
%    Tmap [ 1 x n cell vector of nk x 1 double vectors ]
%    Tlocal [ optional 1 x 4 cell vector of nk x 1 double vectors ]
%       Optional vase Tlocal cell array that's overridden if supplied
% Outputs:
%    Tlocal [ 1 x 4 cell vector of nk x 1 double vectors ]

if nargin < 3
    Tlocalold = [];
end

n = length(Tmap);
sizes = zeros(1,n);
for i = 1:n
    sizes(i) = size(Tmap{i},1);
end

Tlocal = cell(1,n);
for i = 1:n
    
    if isempty(Tlocalold)
        Tlocal_i = zeros(sizes(i),1);
    else
        Tlocal_i = Tlocalold{i};
    end
    
    Tmap_i = Tmap{i};
    for j = 1:sizes(i)
        Tidx = Tmap_i(j);
        if Tidx ~= 0
            Tlocal_i(j) = T(Tidx);
        end
    end
    Tlocal{i} = Tlocal_i;
end

end

function T = mapTlocal2T(Tlocal, Tmap, nT, Told, mode)
% Map local T's in conditions to overall T vector. Allows summing at positions to overlay T.
% Inputs:
%   Tlocal [ 1 x 4 cell vector of nk x 1 double vectors ]
%   Tmap [ 1 x 4 cell vector of nk x 1 double vectors ]
%   nT [ scalar double ]
%       Total number of T parameters
%   Told [ optional nT x 1 double vector ]
%       Optional base vector that's overridden or added to depending
%       on mode
%   mode [ optional string {'override'} ]
%       What newly mapped Tlocal>T should do. If 'override', replace
%       current T with Tlocal value. Used for updating parameters.
%       If 'add', add new value from Tlocal to old T value. Used for
%       gradient calculation.
% Outputs:
%    T [ 1 x nT double vector ]
if nargin < 5
    mode = [];
    if nargin < 4
        Told = [];
    end
end

if isempty(mode)
    mode = 'override';
end

n = length(Tmap);
sizes = zeros(1,n);
for i = 1:n
    sizes(i) = size(Tmap{i},1);
end

if isempty(Told)
    Told = zeros(nT, 1);
end

T = Told;
for i = 1:n
    Tlocal_i = Tlocal{i};
    Tmap_i = Tmap{i};
    
    % Sanity check: Tlocal and Tmap have the same dimensions
    assert(length(Tlocal_i) == length(Tmap_i), 'FitObject:mapTlocal2T: Tlocal and Tmap have different dimensions')
    
    for j = 1:sizes(i)
        Tidx = Tmap_i(j);
        if Tidx ~= 0
            switch mode
                case 'override'
                    T(Tidx) = Tlocal_i(j);
                case 'add'
                    T(Tidx) = T(Tidx) + Tlocal_i(j);
            end
        end
    end
end

end
