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
    % meta properties like partMap
    
   properties
       Name
       modelNames = {}
       conditionNames = {}
       objectiveNames = {}
       nModels = 0
       nConditions = 0;
       nObjectives = 0
       modelVariantsMap = [] % nObjectives double vector of base model index that share param specs
       options
       componentMap = zeros(0,3) % nObjectives x 3 double matrix of [modelIdx, conditionIdx, objectiveIdx]; Note: when parts are turned into objects, get rid of this
       paramsMap = cell(0,4) % nObjectives x 4 cell matrix of [param, seed, input, dose] to T index mapping; 0 indicates not fit; multiple refs to a single T index are summed
       paramsShared = cell(0,4) % % nObjectives x 4 cell matrix of [param, seed, input, dose] shared params index spec
       nT
   end
   
   properties (SetObservable)
       Models
       Conditions
       Objectives
   end
   
   methods
       function this = FitObject(name, models, conditions, objectives, opts)
           % Constructor; Instantiate FitObject.
           
           % Clean up inputs
           if nargin < 5
               opts = [];
               if nargin < 4
                   objectives = [];
                   if nargin < 3
                       conditions = [];
                       if nargin < 2
                           models = [];
                           if nargin < 1;
                               name = [];
                           end
                       end
                   end
               end
           end
           
           % Fill in defaults
           if isempty(name)
               name = 'DefaultFitObject';
           end
           
           % Default options
           opts_.Verbose = 2;
           opts = mergestruct(opts_, opts);
           
           % Input validation and checking
           % TODO
           
           % Attach listeners
           this.addlistener('Models', 'PostSet', @(src,event)this.updateModels(src,event));
           this.addlistener('Conditions', 'PostSet', @(src,event)this.updateConditions(src,event));
           this.addlistener('Objectives', 'PostSet', @(src,event)this.updateObjectives(src,event));
           
           % Assign fields
           this.Name = name;
           this.Models = models;
           this.Conditions = conditions;
           this.Objectives = objectives;
           this.options = opts;
           
       end
       
       function addOption(this, optName, optValue)
           % Add option to fit object
           % 
           % Allowed options:
           %    Aeq
           %    beq
           %
           %TODO: validate inputs
           this.options.(optName) = optValue;
       end
       
       function addModel(this, model, opts)
           if nargin < 3
               opts = [];
           end
           
           nk = model.nk;
           
           % Default options
           opts_.BaseModel = this.nModels + 1; % default is this model is a new type of model
           opts_.LowerBound = 0;
           opts_.UpperBound = Inf;
           opts = mergestruct(opts_, opts);
           
           % Validate options
           LowerBound_ = opts.LowerBound;
           nLowerBounds = length(LowerBound_);
           switch nLowerBounds
               case 1
                   LowerBound_ = repmat(LowerBound_, nk, 1);
               case nk
                   LowerBound_ = vec(LowerBound_);
               otherwise
                   error('FitObject:addModel: Number of lower bounds %d not equal to 1 or number of rate params %d', nLowerBounds, nk)
           end
           opts.LowerBound = LowerBound_;
           
           UpperBound_ = opts.UpperBound;
           nUpperBounds = length(UpperBound_);
           switch nUpperBounds
               case 1
                   UpperBound_ = repmat(UpperBound_, nk, 1);
               case nk
                   UpperBound_ = vec(UpperBound_);
               otherwise
                   error('FitObject:addModel: Number of upper bounds %d not equal to 1 or number of rate params %d', nUpperBounds, nk)
           end
           opts.UpperBound = UpperBound_;
           
           % Annotate models with new information
           model.LowerBound = opts.LowerBound;
           model.UpperBound = opts.UpperBound;
           
           % Add model to fit object
           this.Models = [this.Models; model];
           this.modelVariantsMap = [this.modelVariantsMap; opts.BaseModel];
       end
       
       function addCondition(this, condition, opts)
           % Add experimental condition to fit object.
           % Automatically associates with model with corresponding name (Note: when parts objects, enforce this with typing)
           
           if nargin < 3
               opts = [];
           end
           
           % Check that associated model is already in fit object
           parentModelMask = ismember(this.modelNames, condition.ParentModelName);
           if ~parentModelMask
               error('FitObject:addCondition: Condition parent model %s not present in fit object.', condition.ParentModelName)
           end
           
           % Add associated model index to condition
           condition.ParentModelIdx = find(parentModelMask);
           
           % Add condition to fit object
           this.Conditions = [this.Conditions; condition];
           
       end
       
       function addObjective(this, objective, parentConditions, opts)
           % Add dataset/objective function to fit object.
           % Inputs:
           %    objective:
           %    parentConditions: [ string | cell array of strings {1st condition} ]
           %        names of conditions to associate with
           %    opts:
           
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
           %        .UseSeeds
           %        .UseInputControls
           %        .UseDoseControls
           %
           % By default, the default model, default condition, and all
           % parameters allowed to vary independently option is used.
           % Different models have different param specs (for all k,s,q,h). 
           % Use linear equality constraints to set params shared between models.
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
           
           % Default options
           opts_.UseParams = [];
           opts_.UseSeeds = [];
           opts_.UseInputControls = [];
           opts_.UseDoseControls = [];
           opts = mergestruct(opts_, opts);
           
           % Fix input param specs
           UseParams        = vec(opts.UseParams);
           UseSeeds         = vec(opts.UseSeeds);
           UseInputControls = vec(opts.UseInputControls);
           UseDoseControls  = vec(opts.UseDoseControls);
           
           % Sanity check: make sure useX are valid for model and condition
           assert(this.Models(ismember(this.modelNames, condition.ParentModelName)).nk == length(UseParams), 'FitObject:addFitConditionData: Number of UseParams doesn''t match number of model rate params')
           assert(condition.ns == length(UseSeeds), 'FitObject:addFitConditionData: Number of UseSeeds doesn''t match number of model seed params')
           assert(condition.nq == length(UseInputControls), 'FitObject:addFitConditionData: Number of UseInputControls doesn''t match number of model input control params')
           assert(condition.nh == length(UseDoseControls), 'FitObject:addFitConditionData: Number of UseDoseControls doesn''t match number of model dose control params')
           
           % TODO: API for easy use of linear equality constraints between
           % different models (probably using opts)
           
           % Add UseX for this objective to fit object's collection
           this.paramsShared = [this.paramsShared; {UseParams, UseSeeds, UseInputControls, UseDoseControls}];
           
           % Add dummy model if k is to be fit independently
           modelIdx = this.addModelOnUniqueParam(UseParams, condition.ParentModelName);
           condition.ParentModelName = this.modelNames(modelIdx);
           
           % Update mapping from individual objectives to overall T
           [this.paramsMap, this.nT] = updateParamsMap(this.paramsShared);
           
           % Add components to fit
           this.addCondition(condition, opts);
           this.addObjective(objective, condition.Name, opts);
           
           % Helper functions
           function [paramsMap, nT] = updateParamsMap(paramsShared)
               % Update mapping from individual objectives to overall T
               Toffset = 0;
               
               kShared = paramsShared(:,1);
               kMap = multiModelParamsShared2Map(kShared);
               [kMap, Toffset] = addOffsetToNonzeroCellArray(kMap, Toffset);
               
               sShared = paramsShared(:,2);
               sMap = multiModelParamsShared2Map(sShared);
               [sMap, Toffset] = addOffsetToNonzeroCellArray(sMap, Toffset);
               
               qShared = paramsShared(:,3);
               qMap = multiModelParamsShared2Map(qShared);
               [qMap, Toffset] = addOffsetToNonzeroCellArray(qMap, Toffset);
               
               hShared = paramsShared(:,4);
               hMap = multiModelParamsShared2Map(hShared);
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
           
           function paramsMap = multiModelParamsShared2Map(paramsShared)
               % Inputs:
               %    paramsShared: cell vector nObjective x 1
               % Outputs:
               %    paramsMap: cell vector nObjective x 1
               %
               % Note: logically, rate params are treated differently in
               % different models while seeds, input control params, and dose
               % control params may be more conveniently treated the same in the
               % same conditions, regardless of dummy underlying model for
               % changed k
               
               % Treat different base models separately
               modelTypes = unique(this.modelVariantsMap);
               Toffset = 0;
               paramsMap = cell(size(paramsShared));
               for modelType = modelTypes
                   
                   paramsSharedModel = paramsShared(this.modelVariantsMap == modelType);
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
                   
                   paramsMap(this.modelVariantsMap == modelType) = paramsMapModel;
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
           
       end
       
       function modelIdx = addModelOnUniqueParam(this, useParams, modelName)
           % Add a new "dummy" model if useParams specifies unique rate parameters
           % from any existing for the specified modelName
           % Returns index of model that useParams will refer to, either for an
           % existing model if they match or a new one if it needed to be
           % generated
           % 
           % Note: currently used for rate params but may be adapted for others
           
           % Get appropriate model
           modelType = find(ismember(this.modelNames, modelName));
           
           % Get appropriate rate params specs based on same base model
           paramsSharedModel = this.paramsShared(this.modelVariantsMap == modelType, 1);
           
           % Look for exactly matching params specs
           for i = 1:length(paramsSharedModel)
               if all(paramsSharedModel{i} == useParams)
                   modelIdx = i;
                   return
               end
           end
           
           % No matches - make new model
           model_ = this.getModel(modelName);
           model_.Name = [model_.Name '_gen' num2str(this.nModels)];
           this.Models = [this.Models; model_];
           this.modelVariantsMap = [this.modelVariantsMap; modelType];
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
       
       function summary = summarize(this)
           % Return summary struct of fitting object
       end
       
       %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Called inside optimizer's loop
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function [G, D] = computeObjective(this)
           
           if nargout == 1
           G = 0;
           end
           if nargout == 2
               G = 0;
               D = zeros()
           end
       end
       
       
       function updateParams(this)
%            this.traverse(@)
           
       end
       
       function traverse(this)
           
       end
       
       function Tlocal = Tmain2Tlocal(Tmain)
           
       end
       
       function Tmain = Tlocal2Tmain(Tlocal)
           
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
           % Update nObjectives, objectiveNames, componentMap
           nObjectives_ = length(this.Objectives);
           this.nObjectives = nObjectives_;
           
           objectiveNames_ = cell(nObjectives_,1);
           for i = 1:nObjectives_
               objectiveNames_{i} = this.Objectives(i).Name;
           end
           this.objectiveNames = objectiveNames_;
           
           % Update models < conditions < objectives componentMap which represents the
           % tree as a matrix
           componentMap_ = [];
           for objectiveIdx = 1:this.nObjectives
               conditionName = this.Objectives(objectiveIdx).ParentConditionName;
               conditionIdx = find(ismember(this.conditionNames, conditionName));
               modelName = this.Conditions(conditionIdx).ParentModelName;
               modelIdx = find(ismember(this.modelNames, modelName));
               componentMap_ = [componentMap_; [modelIdx, conditionIdx, objectiveIdx]];
           end
           this.componentMap = componentMap_;
           
       end
       
   end
   
   methods (Static)
       function fit = buildFitObject(models, conditions, objectives, options)
           % Factory function to convert legacy format to fit object
           % Inputs:
           %    models [ struct vector nModels ]
           %        Models to fit. Note: ignores all models except 1st one
           %        (since multi-model fitting wasn't previously implemented).
           %    conditions [ struct vector nConditions ]
           %    objectives [ struct matrix nConditions x nObjectives ]
           %    options [ struct ]
           
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
           
       end
   end
   
   
end
