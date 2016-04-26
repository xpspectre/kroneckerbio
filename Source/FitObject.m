classdef FitObject < handle
    % Specifies fiting scheme for FitObjective
    % Models, Conditions, and Objectives (objective functions) are stored in
    % struct arrays. They can be directly indexed by calling code for
    % convenience.
    % The 1st of each component is the default; other components that are added
    % refer to this one if not otherwise specified.
    % Makes sure entire system is consistent with each component added.
    % Contains all the default options required for fitting (previously in FitObjective)
    
    % TODO/Notes
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
    %       polluting them when they're simply used for simulating - I
    %       needed this functionality for adding stuff to models but worked
    %       hacked in extra fields by making an additional m.Update
    %       function that allows variables to be added.
    
   properties
       Name
       Models
       Conditions
       Objectives
       modelNames = {}
       conditionNames = {}
       objectiveNames = {}
       nModels = 0
       nConditions = 0
       nObjectives = 0
       options
       componentMap = zeros(0,3) % nObjectives x 3 double matrix of [modelIdx, conditionIdx, objectiveIdx]; Note: when parts are turned into objects, get rid of this
       modelClassMap = [] % nObjectives x 1 double matrix of model class index each objective belongs to
       paramMapper % ParamMapper object that holds local condition <-> global fit parameter mapping
       nT
       initializedWrappers = false;
       componentWrappers
%        localWrappers
       nw
   end
   
   methods
       function this = FitObject(name, opts)
           % Constructor; Instantiate FitObject.
           % Inputs:
           %    name [ string {'DefaultFit'} ]
           %        Name of fitting scheme
           %    opts [ struct ]
           %        Options struct allowing the following fields:
           %        .Verbose [ integer scalar {tbd} ]
           %            Larger integers show more debugging output.
           %        .ModelsShareParams [ logical scalar {false} ]
           %            Whether different models (classes of models) share the
           %            same specs for shared/unique params. If true, all param
           %            specs for each of k, s, q, h must be unique among all
           %            models. Facilitates easy shared params between different
           %            models. If false, each model class is treated
           %            independently, individual param specs are
           %            unambiguous, and each k, s, q, h for each model class
           %            has a self contained param spec. Note: false assumes a
           %            single model type.
           %        .Normalized [ logical scalar {true} ]
           %           Indicates if the optimization should be done in log parameters
           %           space
           %    	.UseAdjoint [ logical scalar {true} ]
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
           %        .MaxStepSize [ nonegative scalar {1} ]
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
           %        .ComputeSensPlusOne [ true | {false} ]
           %            Whether to compute the (n+1)-th order sensitivity for an
           %            n-th order calculation. Example: NLME requires 1st order
           %            sensitivities (gradient) for the objective function calculation and
           %            2nd order sensitivities (Hessian) for fitting.
           %        .RunParallelExpts [ true | {false} | scalar positive integer ]
           %            Whether to run experiments in parallel. If true,
           %            initializes a parallel pool using max number of threads
           %            available, or if a scalar positive integer is specified,
           %            using that many threads up to the max number of threads
           %            available. Running experiments in parallel is useful
           %            when the fit needs to integrate many
           %            covariates/experimental conditions/patients for each
           %            step of the optimization.
           %        .GlobalOptimization [ logical scalar {false} ]
           %           Use global optimization in addition to fmincon
           %        .GlobalOpts [ options struct scalar {} ]
           %           TODO: API in progress
           
           % Clean up inputs
           if nargin < 2
               opts = [];
               if nargin < 1
                   name = [];
               end
           end
           
           % Fill in defaults
           if isempty(name)
               name = 'DefaultFit';
           end
           
           % Default options
           opts_.Verbose = 2;
           opts_.ModelsShareParams = false;
           
           opts_.Normalized       = true;
           opts_.UseAdjoint       = true;
           
           opts_.ComplexStep      = false;
           
           opts_.TolOptim         = 1e-5;
           opts_.Restart          = 0;
           opts_.RestartJump      = 0.001;
           opts_.TerminalObj      = -inf;
           
           opts_.MaxStepSize      = 1;
           opts_.Algorithm        = 'active-set';
           opts_.MaxIter          = 1000;
           opts_.MaxFunEvals      = 5000;
           
           opts_.ComputeSensPlusOne = false;
           
           opts_.RunParallelExpts = false;
           
           opts_.Method           = 'fmincon';
           
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
           
           % TODO: validate options
           
           % Assign fields
           this.Name = name;
           this.options = opts;
           this.updateParamMapper;
       end
       
       function addOptions(this, opts)
           % Add options to fit object
           % Inputs:
           %    opts [ struct ]
           %        Options struct allowing the fields specified in FitObject's constructor
           % Side Effects:
           %    Merges in supplied options, overriding existing values/defaults
           %    and discarding invalid options.
           %
           % TODO: validate inputs
           this.options = mergestruct(this.options, opts);
       end
       
       function updateParamMapper(this)
           if this.options.ModelsShareParams
               this.paramMapper = ParamMapperMultiModelType(this.Conditions);
           else
               this.paramMapper = ParamMapperOneModelType(this.Conditions);
           end
       end
       
       %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Main functions for adding fit object components
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function addModel(this, model, opts)
           % Add model to fit object.
           % Inputs:
           %    model [ model struct ]
           %    opts [ struct ]
           %        Options struct allowing the following fields:
           %        .BaseModel [ string ]
           %            Name of model already existing in FitObject. Used
           %            to indicate this model was automatically generated
           %            to serve as a dummy to allow parameters to vary
           %            between conditions based off the same original
           %            model.
           % Side Effects:
           %    Adds model to fit object; updates components map
           
           if nargin < 3
               opts = [];
           end
           
           % Default options
           opts_.BaseModel = model.Name; % default is this model is a new type of model
           opts = mergestruct(opts_, opts);
           
           % Validate and annotate models with new information
           model = model.UpdateExtra(struct('BaseModel', opts.BaseModel));
           
           % Add model to fit object
           this.Models = [this.Models; model];
           
           % Update model components
           this.updateModels;
       end
       
       function addCondition(this, condition, opts)
           % Add experimental condition to fit object. Model condition is based
           % on must already be present, or it throws an error.
           %
           % Inputs:
           %    condition [ condition struct ]
           %    opts [ options struct ]
           %        Options struct allowing the following fields:
           %        .UseParams [ nk x 1 double vector {ones(nk,1)} ]
           %            Specifies rate parameters to be fit, with nonzero
           %            entries indicating fit, with unique/shared between
           %            conditions denoted by different/same values in each row.
           %            Default is all rate params fit, shared between
           %            conditions.
           %        .UseSeeds [ ns x 1 double vector {zeros{ns,1}} ]
           %            Specifies seed parameters to be fit, with nonzero
           %            entries indicating fit, with unique/shared between
           %            conditions denoted by different/same values in each row.
           %            Default is no seed params fit.
           %        .UseInputControls [ nq x 1 double vector {zeros(nq,1)} ]
           %            Specifies input control parameters to be fit, with nonzero
           %            entries indicating fit, with unique/shared between
           %            conditions denoted by different/same values in each row.
           %            Default is no input control params fit.
           %        .UseDoseControls [ nh x 1 double vector {zeros(nh,1)} ]
           %            Specifies dose control parameters to be fit, with nonzero
           %            entries indicating fit, with unique/shared between
           %            conditions denoted by different/same values in each row.
           %            Default is no dose control params fit.
           %        .StartingParams [ nk x 1 double vector {m.k from parent model} ]
           %            Starting rate parameter values for the fit. Overrides
           %            values in the parent model.
           %        .StartingSeeds [ ns x 1 double vector {condition.s} ]
           %            Starting seed parameter values for the fit. Overrides
           %            values in the supplied condition.
           %        .StartingInputControls [ nq x 1 double vector {condition.q} ]
           %            Starting input control parameter values for the fit. Overrides
           %            values in the supplied condition.
           %        .StartingDoseControls [ nh x 1 double vector {condition.h} ]
           %            Starting dose control parameter values for the fit. Overrides
           %            values in the supplied condition.
           %        .ParamLowerBound [ nk x 1 double vector {zeros(nk,1)} ]
           %            Minimum values rate parameters can take in constrained
           %            fit.
           %        .ParamUpperBound [ nk x 1 double vector {inf(nk,1)} ]
           %            Maximum values rate parameters can take in constrained
           %            fit.
           %        .SeedLowerBound [ ns x 1 double vector {zeros(ns,1)} ]
           %        .SeedUpperBound [ ns x 1 double vector {inf(ns,1)} ]
           %        .InputControlLowerBound [ nq x 1 double vector {zeros(nq,1)} ]
           %        .InputControlUpperBound [ nq x 1 double vector {inf(nq,1)} ]
           %        .DoseControlLowerBound [ nh x 1 double vector {zeros(nh,1)} ]
           %        .DoseControlUpperBound [ nh x 1 double vector {inf(nh,1)} ]
           %        .RelTol [ positive scalar double {1e-6} ]
           %            Relative tolerance of the integration
           %        .AbsTol [ positive double vector {1e-9} ]
           %            Absolute tolerance of the integration.
           %            TODO: allow specification of AbsTol for all components
           %            of integration separately, with spec for different
           %            terms in regular integration, forward, and adjoint
           %            sensitivity calculations.
           %        .ParamSpec [ [Name] x [Type,LB,UB,Use] table {empty} ]
           %            Table specifying all parameter bounds and fit specs. The
           %            parameter name is the row and Type, LB, UB, and Use are
           %            columns. Alternative way of specifying fit spec by name
           %            in a nice way.
           %
           % Side Effects:
           %    Adds condition to fit object; updates component map. Annotates
           %    condition.Extra with:
           %        .ParentModelIdx [ scalar double ]
           %            Index of model this condition refers to in fit object
           %        .ParamsSpec [1 x 4 cell vector of nx x 1 double vectors ]
           %            Specifies fit parameters, with unique/shared spec.
           %            Simply a copy of Use* for k, s, q, and h in order in a
           %            cell array.
           %        .Bounds [ 2 x 4 cell vector of nx x 1 double vectors { bounds [0,Inf] for all } ]
           %            Lower (top 1 x 4 cell array) and upper (bottom 1 x 4 cell array)
           %            allowed values for parameters k, s, q, and h. Bounds for
           %            params not fit are ignored.
           %        .RelTol [ positive scalar double {1e-6} ]
           %            Copied from input
           %        .AbsTol [ positive scalar double {1e-9} ]
           %            Copied from input
           %        
           % Note: Automatically associates with model with corresponding
           %    name (when turned into objects, enforce this more strictly)
           
           if nargin < 3
               opts = [];
           end
           
           % Check that associated model is already in fit object
           parentModelName = condition.Extra.ParentModelName;
           parentModelMask = ismember(this.modelNames, parentModelName);
           if ~parentModelMask
               error('FitObject:addCondition: Condition parent model %s not present in fit object.', parentModelName)
           end
           
           % Add associated model index to condition for convenience
           parentModelIdx = find(parentModelMask);
           condition = condition.UpdateExtra(struct('ParentModelIdx', parentModelIdx));
           
           % Correct count of each param type
           nk = this.Models(parentModelIdx).nk;
           ns = condition.ns;
           nq = condition.nq;
           nh = condition.nh;
           
           % Existing params
           kStart = this.Models(parentModelIdx).k;
           sStart = condition.s;
           qStart = condition.q;
           hStart = condition.h;
           
           % Default options
           opts_ = [];
           opts_.UseParams        = ones(nk, 1); % fit all rate params, same in model class
           opts_.UseSeeds         = zeros(ns, 1); % don't fit seeds
           opts_.UseInputControls = zeros(nq, 1); % don't fit input control params
           opts_.UseDoseControls  = zeros(nh, 1); % don't fit dose control params
           
           opts_.StartingParams        = [];
           opts_.StartingSeeds         = [];
           opts_.StartingInputControls = [];
           opts_.StartingDoseControls  = [];
           
           opts_.ParamLowerBound        = zeros(nk,1);
           opts_.ParamUpperBound        = inf(nk,1);
           opts_.SeedLowerBound         = zeros(ns,1);
           opts_.SeedUpperBound         = inf(ns,1);
           opts_.InputControlLowerBound = zeros(nq,1);
           opts_.InputControlUpperBound = inf(nq,1);
           opts_.DoseControlLowerBound  = zeros(nh,1);
           opts_.DoseControlUpperBound  = inf(nh,1);
           
           opts_.RelTol = 1e-6;
           opts_.AbsTol = 1e-9;
           
           opts_.ParamSpec = [];
           
           opts = mergestruct(opts_, opts);
           
           % Parameter names for convenience
           kNames = {this.Models(parentModelIdx).Parameters.Name}';
           sNames = {this.Models(parentModelIdx).Seeds.Name}';
           % TODO: names for input and dose control params?
           
           % Convert opts.ParamSpec and populate corresponding fields
           if ~isempty(opts.ParamSpec)
               % Get param name lists to get indices from
               ps = opts.ParamSpec;
               np = height(ps);
               names = ps.Properties.RowNames;
               for i = 1:np
                   p = ps(i,:);
                   
                   name = names{i};
                   type = p.Type{1};
                   lb   = p.LB;
                   ub   = p.UB;
                   use  = p.Use;
                   
                   switch type
                       case 'k'
                           ind = find(ismember(kNames, name));
                           opts.ParamLowerBound(ind) = lb;
                           opts.ParamUpperBound(ind) = ub;
                           opts.UseParams(ind) = use;
                       case 's'
                           ind = find(ismember(sNames, name));
                           opts.SeedLowerBound(ind) = lb;
                           opts.SeedUpperBound(ind) = ub;
                           opts.UseSeeds(ind) = use;
                       case 'q'
                           error('Not implemented yet, or qNames not specified anywhere')
                       case 'h'
                           error('Not implemented yet, or hNames not specified anywhere')
                       otherwise
                           warning('FitObject:processParamSpec:InvalidParamType', 'Param type %s in opts.ParamSpec not recognized. Skipping.', type)
                   end
               end
           end
           
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
           update = [];
           update.Bounds = bounds;
           condition = condition.UpdateExtra(update);
           
           % Sanity check: make sure useX are valid for model and condition
           assert(nk == length(opts.UseParams), 'FitObject:addCondition: Number of UseParams doesn''t match number of model rate params')
           assert(ns == length(opts.UseSeeds), 'FitObject:addCondition: Number of UseSeeds doesn''t match number of model seed params')
           assert(nq == length(opts.UseInputControls), 'FitObject:addCondition: Number of UseInputControls doesn''t match number of model input control params')
           assert(nh == length(opts.UseDoseControls), 'FitObject:addCondition: Number of UseDoseControls doesn''t match number of model dose control params')
           
           % Validate params specs and add to condition
           update = []; % Need to assemble struct literal because can't call "struct()" fun on cell array
           update.ParamsSpec = {fixParamsSpec(opts.UseParams, nk), ...
                                fixParamsSpec(opts.UseSeeds, ns), ...
                                fixParamsSpec(opts.UseInputControls, nq), ...
                                fixParamsSpec(opts.UseDoseControls, nh)};
           condition = condition.UpdateExtra(update);
           
           % Override starting params that are fit with user specified values
           % If Starting* aren't supplied, default values from underlying model
           % or supplied condition are used.
           kUse = logical(opts.UseParams);
           sUse = logical(opts.UseSeeds);
           qUse = logical(opts.UseInputControls);
           hUse = logical(opts.UseDoseControls);
           
           if ~isempty(opts.StartingParams) && any(kUse)
               kStart(kUse) = validateBounds(opts.StartingParams(kUse), bounds{1,1}(kUse), bounds{2,1}(kUse), kNames(kUse));
           end
           if ~isempty(opts.StartingSeeds) && any(sUse)
               sStart(sUse) = validateBounds(opts.StartingSeeds(sUse), bounds{1,2}(sUse), bounds{2,2}(sUse), sNames(sUse));
           end
           if ~isempty(opts.StartingInputControls) && any(qUse)
               qStart(qUse) = validateBounds(opts.UseInputControls(qUse), bounds{1,3}(qUse), bounds{2,3}(qUse));
           end
           if ~isempty(opts.StartingDoseControls) && any(hUse)
               hStart(hUse) = validateBounds(opts.StartingDoseControls(hUse), bounds{1,4}(hUse), bounds{2,4}(hUse));
           end
           
           if ~isempty(opts.StartingParams) && any(kUse)
               this.Models(parentModelIdx) = this.Models(parentModelIdx).Update(kStart);
           end
           if (~isempty(opts.StartingSeeds) && any(sUse)) || (~isempty(opts.UseInputControls) && any(qUse)) || (~isempty(opts.StartingDoseControls) && any(hUse))
               condition = condition.Update(sStart, qStart, hStart);
           end
           
           % Add tolerances
           condition = condition.UpdateExtra(struct('RelTol', opts.RelTol, ...
               'AbsTol', opts.AbsTol));
           
           % Add condition to fit object
           this.Conditions = [this.Conditions; condition];
           
           % Update conditions maps
           this.updateConditions;
           
       end
       
       function addObjective(this, objective, parentConditions, opts)
           % Add dataset/objective function to fit object.
           %
           % Inputs:
           %    objective [ objective struct ]
           %    parentConditions [ string | cell array of strings {1st condition} ]
           %        names of conditions to associate with
           %    opts [ struct {[]} ]
           %        Options struct allowing the following fields:
           %        .Weight [ double {1} ]
           %            Relative weight of this objective function when multiple
           %            are summed together in fit object. Simply acts as a
           %            multiplier (i.e., no normalization is needed or added)
           %
           % Side Effects:
           %    Adds objective to fit object; updates component maps to be
           %    consistent
           
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
               % Annotate objective with weight
               objective.Weight = opts.Weight;
               
               % Annotate objective with parent condition it's applied to
               %    Objective functions aren't closures so this just works
               parentConditionName = this.conditionNames{conditionsIdxs(i)};
               objective.ParentConditionName = parentConditionName;
               
               % Sanity check: make sure output indices in objective map to actual
               % outputs in model (referenced thru condition)
               parentConditionMask = ismember(this.conditionNames, parentConditionName);
               assert(sum(parentConditionMask) < 2, 'FitObject:addObjective:ConditionNameCollision', 'Multiple conditions have the same name. Make sure conditions have unique names. But we should change this so they''re automatically renamed.')
               parentModelMask = ismember(this.modelNames, this.Conditions(parentConditionMask).Extra.ParentModelName);
               parentModelName = this.Models(parentModelMask).Name;
               outputNames = {this.Models(parentModelMask).Outputs.Name};
               if isnumeric(objective.Outputs)
                   try
                       fitOutputs = outputNames(objective.Outputs);
                   catch
                       error('FitObject:addObjective: Obj fun has outputs not in model')
                   end
               elseif iscell(objective.Outputs)
                   assert(all(ismember(objective.Outputs, outputNames)), 'FitObject:addObjective: Obj fun has outputs not in model');
                   fitOutputs = objective.Outputs;
               elseif ischar(objective.Outputs)
                   assert(ismember(objective.Outputs, outputNames), 'FitObject:addObjective: Obj fun has outputs not in model')
                   fitOutputs = {objective.Outputs};
               end
               
               % Add objective to fit object
               this.Objectives = [this.Objectives; objective];
               
               % If verbose, print brief description of outputs being fit.
               if this.options.Verbose > 2
                   fprintf('Added objective %s that fits outputs %s in condition %s in model %s\n', objective.Name, cellstr2str(fitOutputs), parentConditionName, parentModelName)
               end
           end
           
           % Update objective components
           this.updateObjectives;
           
       end
       
       %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Convenience objective/condition adder for mixed effects participant
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function addFitConditionData(this, objective, condition, opts)
           % Higher level/convenience function that adds a "participant" to fit. Abstracts
           % away some of the built in structure that kroneckerbio uses.
           %
           % Inputs:
           %    objective [ objective struct ]
           %        Dataset/objective function of participant
           %    condition [ condition struct {1st condition} ]
           %        Experimental condition, with seeds, inputs, and doses. 
           %        Selects 1st condition (and it's corresponding model) if not
           %        specified.
           % 	opts [ struct {[]} ]
           %        Options struct allowing the fields specified in
           %        addCondition. Of special note:
           %        .UseParams [ nk x 1 double vector ]
           %            Specifying rate parameters to vary from other conditions
           %            that use the same model type generates a dummy model
           %            with the varying rates.
           %        .ParamSpec
           %
           % Side Effects:
           %    Adds a dummy model if params are allowed to vary between
           %    conditions.
           
           if nargin < 4
               opts = [];
               if nargin < 3
                   condition = [];
               end
           end
           
           if isempty(condition)
               condition = this.Conditions(1);
           end
           
           % Check existing models for existence of valid parent model
           parentModelName = condition.Extra.ParentModelName;
           assert(any(ismember(this.modelNames, parentModelName)), 'FitObject:addFitConditionData: Parent model %s does not exist.', parentModelName)
           
           % If ParamSpec is present, convert it to UseParams for next step
           if isfield(opts, 'ParamSpec')
               kNames = {this.Models(ismember(this.modelNames, parentModelName)).Parameters.Name};
               opts.UseParams = zeros(length(kNames),1);
               % Get param name lists to get indices from
               ps = opts.ParamSpec;
               np = height(ps);
               names = ps.Properties.RowNames;
               for i = 1:np
                   p = ps(i,:);
                   
                   name = names{i};
                   type = p.Type{1};
                   use  = p.Use;
                   
                   if strcmp(type, 'k')
                       opts.UseParams(ismember(kNames, name)) = use;
                   end
                   
               end
           end
           
           % Check if specified fit params match an existing model/condition; if
           %    condition specifies different fit params for an existing model,
           %    make a dummy model. Needed because models specify k inside them
           % If UseParams is not specified, the default is for all params
           %    to be fit with the same vals for each model (and no dummy model is needed) 
           if isfield(opts, 'UseParams')
               
               opts.UseParams = vec(opts.UseParams);
               
               assert(this.Models(ismember(this.modelNames, parentModelName)).nk == length(opts.UseParams), 'FitObject:addFitConditionData: Number of UseParams doesn''t match number of model rate params')
               
               % See if a model with exactly the same params to fit already exists
               %    parentModelIdx = 0 for no valid model; positive integer
               %    for valid model to attach to (usually the 1st
               %    condition)
               parentModelIdx = this.paramMapper.isSharedParamSpec(opts.UseParams);
               
               % See if valid parent model w/o an attached condition exists
               if parentModelIdx == 0 && ~this.hasAttachedCondition(parentModelName)
                   parentModelIdx = find(ismember(this.modelNames, parentModelName));
               end
               
               % Add dummy model if no valid model to attach to (existing parent model attached to other condition)
               %    A model with the m.Extra.DummyName field set uses the
               %    DummyName instead of its regular name when being
               %    referenced by conditions
               if parentModelIdx == 0
                   model_ = this.getModel(parentModelName);
                   dummyName = [model_.Name '_der' num2str(this.nModels)]; % hope this is unique
                   model_ = model_.UpdateExtra(struct('DummyName', dummyName));
                   this.addModel(model_, struct('BaseModel', parentModelName))
                   parentModelName = dummyName;
               end
               
               % Attach condition to model
               condition = condition.UpdateExtra(struct('ParentModelName', parentModelName));
           end
           
           % Add components to fit
           this.addCondition(condition, opts);
           this.addObjective(objective, condition.Name, opts);
       end
       
       function attached = hasAttachedCondition(this, modelName)
           % Get whether or not model specified by modelName has an attached condition
           % Inputs:
           %    modelName [ string ]
           %        Name of model in FitObject to test
           % Outputs:
           %    attached [ scalar boolean ]
           %        Whether or not model has an attached condition
           attached = false;
           if isempty(this.Conditions) % no conditions attached to anything yet
               return
           end
           for i = 1:this.nConditions
               if strcmp(this.Conditions(i).Extra.ParentModelName, modelName)
                   attached = true;
                   return
               end
           end
       end
       
       function model = getModel(this, name)
           % Return model by name.
           % Inputs:
           %    name [ string ]
           %        Name of model to return
           % Outputs:
           %    model [ model struct ]
           %        Model struct with input name
           modelMask = ismember(this.modelNames, name);
           if ~any(modelMask)
               error('FitObject:getModel: Model %s not found', name)
           end
           model = this.Models(modelMask);
       end
       
       function condition = getCondition(this, name)
           % Return condition by name.
           % Inputs:
           %    name [ string ]
           %        Name of condition to return
           % Outputs:
           %    model [ condition struct ]
           %        Condition struct with input name
           conditionMask = ismember(this.conditionNames, name);
           if ~any(conditionMask)
               error('FitObject:getCondition: Condition %s not found', name)
           end
           condition = this.Conditions(conditionMask);
       end
       
       %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Extra optimizer initialization for parallel expt simulation
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function initializeParallelExpts(this)
%            usePar = this.options.RunParallelExpts;
%            if ~usePar % nonzero integers are true
%                warning('FitObject:initializeParallelExpts:NotParallel', 'Running expts in parallel not speified. Ignoring.')
%                return
%            end
%            
%            % Set up parallel pool
%            maxCores = feature('numcores');
%            if islogical(usePar)
%                this.nw = maxCores;
%            elseif isnumeric(usePar)
%                this.nw = min(usePar, maxCores);
%            else
%                error('FitObject:initializeParallelExpts:InvalidOpt', 'RunParallelExpts must be true or a number of threads.')
%            end
%            pool = gcp('nocreate'); % if no pool, do not create new one.
%            if isempty(pool)
%                parpool('local', this.nw);
%            end
           
           componentMaps = splitComponentMap(this.componentMap, this.nw);
%            this.componentWrappers = cell(this.nw,1);
%            for i = 1:this.nw
%                this.componentWrappers{i} = WorkerObjWrapper(@getComponentChunks, {componentMaps{i}, this.Models, this.Conditions, this.Objectives});
%            end
           % Debugging: view the wrappers locally
%            this.localWrappers = cell(this.nw,1);
%            for i = 1:this.nw
%                this.localWrappers{i} = getComponentChunks(componentMaps{i}, this.Models, this.Conditions, this.Objectives);
%            end
            
           
       end
       
       function computeParallel(this, fun)
           % Contains the only SPMD blocks
           
           % Initialize if not done yet
           if ~this.initializedWrappers
               usePar = this.options.RunParallelExpts;
               if ~usePar % nonzero integers are true
                   warning('FitObject:initializeParallelExpts:NotParallel', 'Running expts in parallel not speified. Ignoring.')
                   return
               end
               
               % Set up parallel pool
               maxCores = feature('numcores');
               if islogical(usePar)
                   this.nw = maxCores;
               elseif isnumeric(usePar)
                   this.nw = min(usePar, maxCores);
               else
                   error('FitObject:initializeParallelExpts:InvalidOpt', 'RunParallelExpts must be true or a number of threads.')
               end
               pool = gcp('nocreate'); % if no pool, do not create new one.
               if isempty(pool)
                   parpool('local', this.nw);
               end
               
               componentMaps = splitComponentMap(this.componentMap, this.nw);
               spmd
                   components = getComponentChunks(componentMaps{labindex}, this.Models, this.Conditions, this.Objectives);
               end
               
               this.initializedWrappers = true;
           end
           
           spmd
               components.nConditions
               components.componentMap
               
               Gs = fun(components, this.options)
           end
           
       end
       
       %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Called inside optimizer's loop
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function [G, D, H] = computeObjective(this, conditions)
           % Compute combined objective function value, gradient, and
           % Hessian of this FitObject.
           %
           % Inputs:
           %    conditions [ double vector {1:this.nConditions} ]
           %        Indices of conditions to compute quantities for.
           %        Default is to compute using all conditions. When using
           %        RunParallelExpts, all conditions are always calculated.
           %        TODO: make this work for parallel expts too.
           %
           % Outputs:
           %    G [ double scalar ]
           %        Objective function value from sum of all objectives in fit
           %        object with current parameters
           %    D [ nT x 1 double vector ]
           %        Gradient w.r.t. all parameters in Theta for fmincon
           %    H [ nT x nT double matrix ]
           %        Hessian w.r.t. all parameters in Theta for fmincon
           %
           % Note: If this.options.ComputeSensPlusOne = true, Hessian w.r.t.
           %    parameters to be optimized on cannot be computed.
           % TODO: Ability to use finite difference approximations if specified
           
           if nargin < 2
               conditions = [];
           end
           if isempty(conditions)
               conditions = 1:this.nConditions;
           end
           nCon = length(conditions);
           
           if this.options.ComputeSensPlusOne
               assert(nargout < 3, 'KroneckerBio:FitObject:computeObjective:HessianOrderTooHigh', 'Hessian cannot be computed since ComputeSensPlusOne requires 3rd order sensitivities.')
           end
           
           if nargout == 1 % objective function value
               
               if this.options.RunParallelExpts;
                   if ~this.initializedWrappers; this.initializeParallelExpts; this.initializedWrappers = true; end
                   spmd % Run the experiments in each chunk
                       components = this.componentWrappers{labindex}.Value;
                       nCon = components.nConditions;
                       
                       Gs_ = zeros(nCon,1);
                       for i = 1:nCon
                           modelIdx = components.Conditions(i).Extra.ParentModelIdx;
                           objectives = components.Objectives(components.componentMap(:,2) == i)';
                           
                           if this.options.ComputeSensPlusOne
                               Gs_(i) = computeObjGrad(components.Models(modelIdx), components.Conditions(i), objectives, this.options);
                           else
                               Gs_(i) = computeObj(components.Models(modelIdx), components.Conditions(i), objectives, this.options);
                           end
                       end
                   end
                   Gs = joinList(gather(Gs_));
               else
                   Gs = zeros(nCon,1);
                   for i = 1:nCon
                       iCon = conditions(i);
                       modelIdx = this.Conditions(iCon).Extra.ParentModelIdx;
                       objectives = this.Objectives(this.componentMap(:,2) == iCon)';
                       
                       if this.options.ComputeSensPlusOne
                           Gs(i) = computeObjGrad(this.Models(modelIdx), this.Conditions(iCon), objectives, this.options);
                       else
                           Gs(i) = computeObj(this.Models(modelIdx), this.Conditions(iCon), objectives, this.options);
                       end
                   end
               end
               G = sum(Gs);
               
           end
           
           if nargout == 2 % gradient
               
               if this.options.RunParallelExpts;
                   if ~this.initializedWrappers; this.initializeParallelExpts; this.initializedWrappers = true; end
                   spmd % Run the experiments in each chunk
                       components = this.componentWrappers{labindex}.Value;
                       nCon = components.nConditions;
                       
                       Gs_ = zeros(nCon,1);
                       Ds_ = cell(nCon,1);
                       for i = 1:nCon
                           modelIdx = components.Conditions(i).Extra.ParentModelIdx;
                           objectives = components.Objectives(components.componentMap(:,2) == i)';
                           
                           if this.options.ComputeSensPlusOne
                               [Gs_(i), Ds_(i)] = computeObjHess(components.Models(modelIdx), components.Conditions(i), objectives, this.options);
                           else
                               [Gs_(i), Ds_(i)] = computeObjGrad(components.Models(modelIdx), components.Conditions(i), objectives, this.options);
                           end
                       end
                   end
                   Gs = joinList(gather(Gs_));
                   Ds = joinList(gather(Ds_));
               else
                   Gs = zeros(nCon,1);
                   Ds = cell(nCon,1);
                   for i = 1:nCon
                       iCon = conditions(i);
                       modelIdx = this.Conditions(iCon).Extra.ParentModelIdx;
                       objectives = this.Objectives(this.componentMap(:,2) == iCon)';
                       
                       if this.options.ComputeSensPlusOne
                           [Gs(i), Ds(i)] = computeObjHess(this.Models(modelIdx), this.Conditions(iCon), objectives, this.options);
                       else
                           [Gs(i), Ds(i)] = computeObjGrad(this.Models(modelIdx), this.Conditions(iCon), objectives, this.options);
                       end
                   end
               end
               G = sum(Gs);
               D = this.paramMapper.Tlocal2T(Ds, 'sum');
               
           end
           
           if nargout == 3 % Hessian, normally not called during fitting
               
               if this.options.RunParallelExpts;
                   if ~this.initializedWrappers; this.initializeParallelExpts; this.initializedWrappers = true; end
                   spmd % Run the experiments in each chunk
                       components = this.componentWrappers{labindex}.Value;
                       nCon = components.nConditions;
                       
                       Gs_ = zeros(nCon,1);
                       Ds_ = cell(nCon,1);
                       Hs_ = cell(nCon,1);
                       for i = 1:nCon
                           modelIdx = components.Conditions(i).Extra.ParentModelIdx;
                           objectives = components.Objectives(components.componentMap(:,2) == i)';
                           
                           if this.options.ComputeSensPlusOne
                               error('FitObject:computeObjective:tooHighSensOrder', 'ComputeSensPlusOne = true was specified for an obective Hessian calculation and 3rd order sensitivities are not implemented.')
                           else
                               [Gs_(i), Ds_(i), Hs_(i)] = computeObjHess(components.Models(modelIdx), components.Conditions(i), objectives, this.options);
                           end
                       end
                   end
                   Gs = joinList(gather(Gs_));
                   Ds = joinList(gather(Ds_));
                   Hs = joinList(gather(Hs_));
               else
                   Gs = zeros(nCon,1);
                   Ds = cell(nCon,1);
                   Hs = cell(nCon,1);
                   for i = 1:nCon
                       iCon = conditions(i);
                       modelIdx = this.Conditions(iCon).Extra.ParentModelIdx;
                       objectives = this.Objectives(this.componentMap(:,2) == iCon)';
                       
                       if this.options.ComputeSensPlusOne
                           error('FitObject:computeObjective:tooHighSensOrder', 'ComputeSensPlusOne = true was specified for an obective Hessian calculation and 3rd order sensitivities are not implemented.')
                       else
                           [Gs(i), Ds(i), Hs(i)] = computeObjHess(this.Models(modelIdx), this.Conditions(iCon), objectives, this.options);
                       end
                   end
               end
               G = sum(Gs);
               D = this.paramMapper.Tlocal2T(Ds, 'sum');
               H = this.paramMapper.Tlocal2T(Hs, 'sum');
               
           end
       end
       
       function [F, All] = computeInformation(this, dxdTSol)
           % Compute Fisher information matrix of this FitObject.
           % Optionally accepts sensitivities in case you don't wan't to
           % recalculate them and optionally returns individual FIMs for
           % each condition.
           %
           % Inputs:
           %    dxdTSol [ nCon x 1 cell array of sensitivity solution structs {} ]
           %        Optional sensitivity solution struct(s) from
           %        integrateAllSens if you don't want to recalculate them
           %        in this function. The number and ordering of structs
           %        must match those in this FitObject. 
           %
           % Outputs:
           %   F [ double matrix ]
           %       Combined fisher information matrix. The FIM is the sum of all
           %       fisher information matrices assuming there is no
           %       covariance between errors in seperate experiments.
           %   All [ nCon x 1 cell array of double matrices ]
           %       The individual FIMs for each experiment
           
           if nargin < 2
               dxdTSol = [];
           end
           
           nCon = this.nConditions;
           if isempty(dxdTSol)
               dxdTSol = cell(nCon,1);
           end
           
           if this.options.RunParallelExpts;
               if ~this.initializedWrappers; this.initializeParallelExpts; this.initializedWrappers = true; end
               dxdTSolChunks = splitList(dxdTSol, this.nw); % Divide up dxdTSol if provided, cell array of cell arrays
               
               spmd % Run the experiments in each chunk
                   components = this.componentWrappers{labindex}.Value;
                   nCon = components.nConditions;
                   dxdTSols = dxdTSolChunks{labindex};
                   
                   Fs_ = cell(nCon,1);
                   for i = 1:nCon
                       modelIdx = components.Conditions(i).Extra.ParentModelIdx;
                       objectives = components.Objectives(components.componentMap(:,2) == i)';
                       
                       Fs_(i) = computeObjInfo(components.Models(modelIdx), components.Conditions(i), objectives, this.options, dxdTSols{i});
                   end
               end
               Fs = joinList(gather(Fs_));
           else
               
               Fs = cell(nCon,1);
               for i = 1:nCon
                   modelIdx = this.Conditions(i).Extra.ParentModelIdx;
                   objectives = this.Objectives(this.componentMap(:,2) == i)';
                   
                   Fs(i) = computeObjInfo(this.Models(modelIdx), this.Conditions(i), objectives, this.options, dxdTSol{i});
               end
           end
           F = this.paramMapper.Tlocal2T(Fs, 'sum');
           
           % Return all individual condition FIMs
           if nargout == 2
               nCon = this.nConditions; % get this again if turned into a Composite object inside spmd
               All = cell(nCon,1);
               for i = 1:nCon
                   All{i} = this.paramMapper.Tlocal2T(Fs(i), 'sum');
               end
           end
           
       end
       
       function [T, details] = collectParams(this)
           % Collect all params from fit object into master T. Returns details
           %    on collected parameters as second output arg
           %
           % Outputs:
           %    T [ nT x 1 double vector ]
           %        Vector of main parameters optimized by fmincon.
           %    details [ nT x [Name,Type,Ind] Table ]
           %        Table of parameter name, type {'k','s','q','h'} and the
           %        model (for k) or experiment (for s,q,h) it came from
           
           % Put params in same format as param mapper
           if this.options.RunParallelExpts;
               if ~this.initializedWrappers; this.initializeParallelExpts; this.initializedWrappers = true; end
               spmd
                   components = this.componentWrappers{labindex}.Value;
                   nCon = components.nConditions;
                   Tlocals_ = cell(nCon,1);
                   for i = 1:nCon
                       parentModelIdx = components.Conditions(i).Extra.ParentModelIdx;
                       Tlocals_{i} = {components.Models(parentModelIdx).k; ...
                           components.Conditions(i).s; ...
                           components.Conditions(i).q; ...
                           components.Conditions(i).h};
                   end
               end
               Tlocals = joinList(gather(Tlocals_));
           else
               nCon = this.nConditions;
               Tlocals = cell(nCon,1);
               for i = 1:nCon
                   parentModelIdx = this.Conditions(i).Extra.ParentModelIdx;
                   Tlocals{i} = {this.Models(parentModelIdx).k; ...
                       this.Conditions(i).s; ...
                       this.Conditions(i).q; ...
                       this.Conditions(i).h};
               end
           end
           [T, details_] = this.paramMapper.Tlocal2T(Tlocals);
           
           % Process details_ into nicer details
           if nargout == 2
               paramTypes = {'k','s','q','h'};
               nT = numel(T);
               
               Name = cell(nT,1);
               Type = cell(nT,1);
               Ind = cell(nT,1);
               
               for i = 1:nT
                   Ind{i} = details_(i,3);
                   Type{i} = paramTypes{details_(i,1)};
                   switch Type{i}
                       case 'k'
                           kNames = {this.Models(Ind{i}).Parameters.Name};
                           Name{i} = kNames{details_(i,2)};
                       case 's'
                           parentModelIdx = this.Conditions(Ind{i}).Extra.ParentModelIdx;
                           sNames = {this.Models(parentModelIdx).Seeds.Name};
                           Name{i} = sNames{details_(i,2)};
                       otherwise
                           error('q and h and other param types not implemented yet')
                   end
               end
               
               details = table(Name,Type,Ind);
           end
           
       end
       
       function updateParams(this, T)
           % Update all params in fit object with master T
           % Inputs:
           %    T [ nT x 1 double vector ]
           %        Vector of main parameters optimized by fmincon. Must be
           %        supplied in absolute units (e.g, if options.Normalized, 
           %        exponentiate in calling function).
           % Side Effects:
           %    Updates all models and conditions with T provided
           
           % Get params in param mapper format
           Tlocals = this.paramMapper.T2Tlocal(T);
           
           % Update conditions with Tlocals
           if this.options.RunParallelExpts;
               if ~this.initializedWrappers; this.initializeParallelExpts; this.initializedWrappers = true; end
               TlocalsChunks = splitList(Tlocals, this.nw);
               
               spmd
                   components = this.componentWrappers{labindex}.Value;
                   Tlocals = TlocalsChunks{labindex};
                   
                   for i = 1:components.nConditions
                       parentModelIdx = components.Conditions(i).Extra.ParentModelIdx;
                       
                       % Paste updated param values over old values containing params not fit
                       Tlocalold = {components.Models(parentModelIdx).k, components.Conditions(i).s, components.Conditions(i).q, components.Conditions(i).h};
                       Tlocal = combineCellMats(Tlocalold, Tlocals{i});
                       [k,s,q,h] = deal(Tlocal{:});
                       
                       components.Models(parentModelIdx) = components.Models(parentModelIdx).Update(k);
                       components.Conditions(i) = components.Conditions(i).Update(s, q, h);
                       this.componentWrappers{labindex}.Value = components; %%TODO%% Doesn't work - read only
                   end
               end
           else
               for i = 1:this.nConditions
                   parentModelIdx = this.Conditions(i).Extra.ParentModelIdx;
                   
                   % Paste updated param values over old values containing params not fit
                   Tlocalold = {this.Models(parentModelIdx).k, this.Conditions(i).s, this.Conditions(i).q, this.Conditions(i).h};
                   Tlocal = combineCellMats(Tlocalold, Tlocals{i});
                   [k,s,q,h] = deal(Tlocal{:});
                   
                   this.Models(parentModelIdx) = this.Models(parentModelIdx).Update(k);
                   this.Conditions(i) = this.Conditions(i).Update(s, q, h);
               end
           end
       end
       
       function [LowerBound, UpperBound] = collectBounds(this)
           % Collect all bounds from fit object into master list
           % Used by normalize and for fmincon options
           % Outputs:
           %    LowerBound [ nT x 1 double vector ]
           %        Vector of lower bounds for params optimized by fmincon
           %    UpperBound [ nT x 1 double vector ]
           %        Vector of upper bounds for params optimized by fmincon
           
           % Put bounds in same format as param mapper
           LowerBounds = cell(this.nConditions,1);
           UpperBounds = cell(this.nConditions,1);
           for i = 1:this.nConditions
               LowerBounds{i} = this.Conditions(i).Extra.Bounds(1,:)';
               UpperBounds{i} = this.Conditions(i).Extra.Bounds(2,:)';
           end
           LowerBound = this.paramMapper.Tlocal2T(LowerBounds);
           UpperBound = this.paramMapper.Tlocal2T(UpperBounds);
       end
       
       %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Update fit object components to ensure consistency
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function updateModels(this)
           % Update fit object fields when model is added
           % Uses the model.Extra.DummyName field if present, reflecting an
           %    automatically added model
           this.nModels = length(this.Models);
           this.modelNames = getModelNames(this.Models);
       end
       
       function updateConditions(this)
           % Update fit object fields when condition is added
           this.nConditions = length(this.Conditions);
           this.conditionNames = getConditionNames(this.Conditions);
       end
       
       function updateObjectives(this)
           % Update objectives and all mapping fields for traversing components
           % Called after an objective is added
           nObjectives_ = length(this.Objectives);
           this.nObjectives = nObjectives_;
           
           objectiveNames_ = cell(nObjectives_,1);
           for i = 1:nObjectives_
               objectiveNames_{i} = this.Objectives(i).Name;
           end
           this.objectiveNames = objectiveNames_;
           
           % Update component map
           this.componentMap = buildComponentMap(this.Models, this.Conditions, this.Objectives);
           
           % Update mapping of objectives to their base model classes
           this.updateModelClassMap;
           
           % Update condition <-> fit parameters mapping
           this.updateParamMapper;
       end
       
       function updateModelClassMap(this)
           % Update mapping of objectives to their base model classes
           modelClassMap_ = zeros(this.nObjectives, 1);
           baseModels = cell(0,1);
           for i = 1:this.nObjectives
               baseModel = this.Models(this.componentMap(i,1)).Extra.BaseModel;
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
       
   end
   
   
   %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Static factory functions for convenience/basic backwards compatibility
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   methods (Static)
       function fit = buildFitObject(model, conditions, objectives, opts)
           % Factory function to convert legacy format to fit object. Note that
           % this provides somewhat restricted functionality compared to
           % building the FitObject directly.
           %
           % Inputs:
           %    model [ scalar model struct ]
           %        A single kroneckerbio model
           %    conditions [ nCon x 1 conditions struct vector ]
           %        Experimental conditions under which the model is run
           %    objectives [ nCon x nObj objectives struct matrix ]
           %        Objective functions under which the experimental conditions
           %        are run. Each row of objectives has a corresponding row in
           %        conditions.
           %    opts [ options struct ]
           %        Options struct allowing the following fields:
           %        .UseParams [ nk x 1 logical vector {true(nk,1)} ]
           %            Specifies which rate parameters are fit
           %        .UseSeeds [ ns x 1 logical vector {false(ns,1)} | nCon x 1 cell vector of ns_coni logical vectors]
           %            Specifies which seeds are fit for a single condition or if
           %            all conditions have the same number of seeds. A cell
           %            vector can be specified, for which each element
           %            corresponds to a condition.
           %        .UseInputControls [ nq x 1 logical vector {false(nq,1)} | nCon x 1 cell vector of nq_coni logical vectors ]
           %            Specifies which input control parameters are fit for a single condition or if
           %            all conditions have the same number of seeds. A cell
           %            vector can be specified, for which each element
           %            corresponds to a condition.
           %        .UseDoseControls [ nh x 1 logical vector {false(nh,1)} | nCon x 1 cell vector of nh_coni logical vectors ]
           %            Specifies which dose control parameters are fit for a single condition or if
           %            all conditions have the same number of seeds. A cell
           %            vector can be specified, for which each element
           %            corresponds to a condition.
           %        .ObjWeights [ nCon x nObj real matrix {ones(nCon,nObj)} ]
           %           Applies a post evaluation weight on each objective function
           %           in terms of how much it will contribute to the final objective
           %           function value
           %        .RelTol [ scalar nonnegative double {1e-6} ]
           %            Relative tolerance for integration
           %        .AbsTol [ nCon x 1 cell vector of AbsTol structs | nCon x 1 cell vector of double vectors {1e-9 for everything} ]
           %            Absolute tolerance for integration
           %        .Normalized [ scalar logical {true} ]
           %            Whether to work with parameters in log space
           %        .UseAdjoint [ scalar logical {true} ]
           %            Whether to use adjoint (forward if false) method to
           %            calculate sensitivities
           %        .ComplexStep [ scalar logical {false} ]
           %            Whether to use complex step difference (finite
           %            difference if false) when finite difference
           %            approximation is specified
           %        .Verbose [ scalar nonnegative integer {1} ]
           %            Higher verbosity indicates more printed debug output
           %
           % Outputs:
           %    fit [ scalar FitObject ]
           %        Fitting scheme
           
           if nargin < 4
               opts = [];
           end
           
           % Basic sanity checks
           assert(isscalar(model), 'FitObject:buildFitObject:NonscalarModel', 'The model structure must be scalar')
           nCon = length(conditions);
           [nObjCon, nObj] = size(objectives);
           assert(nObjCon == nCon, 'FitObject:buildFitObject:ObjDimMismatch', 'The number of objective function rows %g doesn''t match the number of experimental conditions %g', nObjCon, nCon)
           
           nk = model.nk;
           ns = model.ns;
           nqs = [conditions.nq]';
           nhs = [conditions.nh]';
           useInputControls = cell(nCon,1);
           useDoseControls = cell(nCon,1);
           for i = 1:nCon
               useInputControls{i} = false(nqs(i),1);
               useDoseControls{i}  = false(nhs(i),1);
           end
           
           % Set default options
           opts_ = [];
           opts_.UseParams        = true(nk,1);
           opts_.UseSeeds         = false(ns,1);
           opts_.UseInputControls = useInputControls;
           opts_.UseDoseControls  = useDoseControls;
           opts_.ObjWeights       = ones(size(objectives));
           opts_.RelTol = 1e-6;
           opts_.AbsTol = 1e-9;
           opts_.Normalized = true;
           opts_.UseAdjoint = true;
           opts_.CompexStep = false;
           opts_.Verbose = 1;
           opts_.continuous = false;
           opts = mergestruct(opts_, opts);
           
           % Convert old format options to new format ones
           useParams = (1:nk)'.*fixUseParams(opts.UseParams, nk);
           useSeeds = repmat((1:ns), nCon, 1)'.*fixUseSeeds(opts.UseSeeds, ns, nCon); % all seeds are unique
           useInputControls = fixUseControls(opts.UseInputControls, nCon, cat(1,conditions.nq));
           useDoseControls = fixUseControls(opts.UseDoseControls, nCon, cat(1,conditions.nh));
           
           % Standardize AbsTol
           absTol = repmat({opts.AbsTol}, nCon, 1);
           
           % Make new fit object
           fit = FitObject('ConvertedFit');
           
           % Add rest of options (invalid options ignored when added to fit object)
           fit.addOptions(opts);
           
           % Add model
           fit.addModel(model);
           
           % Add conditions
           for i = 1:nCon
               conditionOpts = [];
               conditionOpts.UseParams        = useParams;
               conditionOpts.UseSeeds         = useSeeds(:,i);
               conditionOpts.UseInputControls = useInputControls{i};
               conditionOpts.UseDoseControls  = useDoseControls{i};
               conditionOpts.AbsTol           = absTol{i};
               fit.addCondition(conditions(i), conditionOpts);
           end
           
           % Add objectives
           for i = 1:nCon
               conditionName = conditions(i).Name;
               for j = 1:nObj
                   obj = objectives(i,j);
                   if ~strcmpi(obj.Type, 'Objective.Data') % only add non-objectiveZeros, assuming this generic name isn't overridden
                       fit.addObjective(obj, conditionName);
                   end
               end
           end
           

       end
   end
   
   
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

function out = validateBounds(in, lb, ub, names)
% Check bounds on parameters: returning it if it satisfies them, setting equal to
% appropriate bound if it falls outside and throwing warning. Accepts vectors of
% in, lb, and ub. Throws error if input vectors' lengths don't match.
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
assert(nin == nlb, 'FitObject:validateBounds:lbLengthMismatch', 'Fit param vector length %d and lower bounds length %d don''t match', nin, nlb);
assert(nin == nub, 'FitObject:validateBounds:ubLengthMismatch', 'Fit param vector length %d and upper bounds length %d don''t match', nin, nub);
if ~isempty(names)
    nNames = length(names);
    assert(nin == nNames, 'FitObject:validateBounds:nameLengthMismatch', 'Fit param vector length %d and names length %d don''t match', nin, nNames);
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
        warning('FitObject:validateBounds: fit param %s = %d below lower bound = %d. Setting to lower bound', name, in(lInd), lb(lInd))
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

function fixedParamsSpec = fixParamsSpec(ParamsSpec, n)
nParams = length(ParamsSpec);
switch nParams
    case n
        fixedParamsSpec = vec(ParamsSpec);
    otherwise
        error('FitObject:fixParamsSpec: Number of params specs %d not equal to number params %d', nParams, n)
end
end

function modelNames = getModelNames(Models)
% Gets model names, using m.Extra.DummyName if present
nModels = length(Models);
modelNames = cell(nModels,1);
for i = 1:nModels
    if isfield(Models(i).Extra, 'DummyName')
        modelNames{i} = Models(i).Extra.DummyName;
    else
        modelNames{i} = Models(i).Name;
    end
end
end

function parentModelIdx = getParentModelIdx(Condition, Models)
parentModelName = Condition.Extra.ParentModelName;
modelNames = getModelNames(Models);
parentModelIdx = find(ismember(modelNames, parentModelName));
assert(numel(parentModelIdx) == 1, 'FitObject:getParentModelIdx:MultipleParentModelsFound', 'Multiple parent models were found for Condition based on ParentModelName')
end

function conditionNames = getConditionNames(Conditions)
% Gets condition names
nConditions = length(Conditions);
conditionNames = cell(nConditions,1);
for i = 1:nConditions
    conditionNames{i} = Conditions(i).Name;
end
end

function componentMap = buildComponentMap(Models, Conditions, Objectives)
% Build models < conditions < objectives componentMap which represents the
%   tree as a matrix
nObjectives = length(Objectives);
conditionNames = getConditionNames(Conditions);
modelNames = getModelNames(Models);
componentMap = zeros(nObjectives, 3);
for i = 1:nObjectives
    conditionName = Objectives(i).ParentConditionName;
    conditionIdx = find(ismember(conditionNames, conditionName));
    modelName = Conditions(conditionIdx).Extra.ParentModelName;
    modelIdx = find(ismember(modelNames, modelName));
    componentMap(i,:) = [modelIdx, conditionIdx, i];
end
end

function components = getComponentChunks(componentMap, Models, Conditions, Objectives)
% Returns self-contained struct of components for use with WorkerObjWrapper
% COW for Models, Components, and Objectives means this should be fast
components = [];
components.nConditions = size(componentMap,1);
components.Models = Models(componentMap(:,1));
components.Conditions = Conditions(componentMap(:,2));
for i = 1:components.nConditions % Reassign Condition ParentModelIdx to point to the proper Model position
    components.Conditions(i) = components.Conditions(i).UpdateExtra(struct('ParentModelIdx', getParentModelIdx(components.Conditions(i), components.Models)));
end
components.Objectives = Objectives(componentMap(:,3));
components.componentMap = buildComponentMap(components.Models, components.Conditions, components.Objectives);
end
