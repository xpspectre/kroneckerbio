function varargout = computeObjAll(m, con, obj, opts, mode, varargin)
% Compute various objective function evaulations. Expects input args in new
% format.
%
% Usage:
%   G = computeObjAll(m, con, obj, opts, 'val')
%   [G, D] = computeObjAll(m, con, obj, opts, 'grad')
%   [G, D, H] = computeObjAll(m, con, obj, opts, 'hess')
%   F = computeObjAll(m, con, obj, opts, 'info')
%   [F, FAll] = computeObjAll(m, con, obj, opts, 'info')
%   [F, FAll] = computeObjAll(m, con, obj, opts, 'info', dxdTSol)
%   p = computeObjAll(m, con, obj, opts, 'prob')
%   logp = computeObjAll(m, con, obj, opts, 'logp')
%
% Inputs:
%	m [ nModel x 1 model struct vector ]
%       The KroneckerBio model that will be simulated
% 	con [ nCon x 1 experiment struct vector ]
%       The experimental conditions under which the model will be simulated
%	obj [ nObj x 1 objective struct vector ]
%       The objective structures defining the objective functions to be
%       evaluated.
%	opts [ options struct scalar ]
%       Options struct with just global opts or detailed options from
%       BuildFitOpts.
%   mode [ 'val' | 'grad' | 'hess' | 'info' | 'prob' | 'logp' ]
%       The particular objective function evaluation desired.
%   dxdTSol [ nCon x 1 cell array of sensitivity solution structs {} ]
%       Optional sensitivity solution struct(s) from
%       integrateAllSens if you don't want to recalculate them
%       in this function. The number and ordering of structs
%       must match those in con/opts.
%
% Outputs:
%	G [ double scalar ]
%       Objective function value from sum of all objectives in fit
%       object with current parameters
%   D [ nT x 1 double vector ]
%       Gradient w.r.t. all parameters in Theta for fmincon
%   H [ nT x nT double matrix ]
%       Hessian w.r.t. all parameters in Theta for fmincon
%   F [ double matrix ]
%       Combined fisher information matrix. The FIM is the sum of all
%       fisher information matrices assuming there is no
%       covariance between errors in seperate experiments.
%   FAll [ nCon x 1 cell array of double matrices ]
%       The individual FIMs for each experiment
%   p [ real nonnegative scalar ]
%       The product of all objective function probabilities
%   logp [ real nonnegative scalar ]
%       The log of the product of all objective function probabilities
%
% Notes:
% - The log-likelihood logp (the negative of it anyway) can basically serve as
%   the objective function value G for all intents and purposes (and more for
%   true MLE usage)

nCon = length(con);
ComponentMap = opts.fit.ComponentMap;

if opts.ComputeSensPlusOne
    assert(nargout < 3, 'KroneckerBio:computeObjAll:HessianOrderTooHigh', 'Hessian cannot be computed since ComputeSensPlusOne requires 3rd order sensitivities.')
end

if strcmp(mode, 'info')
    if length(varargin) == 1 && ~isempty(varargin{1})
        dxdTSol = varargin{1};
        assert(length(dxdTSol) == nCon);
    else
        dxdTSol = cell(nCon,1);
    end
end

Gs = zeros(nCon,1);
Ds = cell(nCon,1);
Hs = cell(nCon,1);
Fs = cell(nCon,1);
ps = zeros(nCon,1);
logps = zeros(nCon,1);
for iCon = 1:nCon
    m_i = m(ComponentMap{iCon,1});
    con_i = con(ComponentMap{iCon,2});
    obj_i = obj([ComponentMap{iCon,3}])';
    
    opts_i = [];
    opts_i.Verbose = opts.Verbose;
    opts_i.Normalized = opts.Normalized;
    opts_i.UseAdjoint = opts.UseAdjoint;
    opts_i.RelTol = opts.RelTol;
    opts_i.AbsTol = opts.fit.AbsTol{iCon};
    opts_i.paramSpec = opts.fit.ParamSpec{iCon}; % param spec for just this condition
    opts_i.UseParams = logical(opts.fit.UseParams{iCon});
    opts_i.UseSeeds = logical(opts.fit.UseSeeds{iCon});
    opts_i.UseInputControls = logical(opts.fit.UseInputControls{iCon});
    opts_i.UseDoseControls = logical(opts.fit.UseDoseControls{iCon});
    opts_i.Continuous = opts.fit.Continuous(iCon);
    opts_i.ObjWeights = opts.fit.ObjWeights{iCon};
    
    switch mode
        case 'val'
            if opts.ComputeSensPlusOne
                Gs(iCon) = computeObjGrad(m_i, con_i, obj_i, opts_i);
            else
                Gs(iCon) = computeObj(m_i, con_i, obj_i, opts_i);
            end
        case 'grad'
            if opts.ComputeSensPlusOne
                [Gs(iCon), Ds(iCon)] = computeObjHess(m_i, con_i, obj_i, opts_i);
            else
                [Gs(iCon), Ds(iCon)] = computeObjGrad(m_i, con_i, obj_i, opts_i);
            end
        case 'hess'
            if opts.ComputeSensPlusOne
                error('FitObject:computeObjAll:HessianOrderTooHigh', 'ComputeSensPlusOne = true was specified for an obective Hessian calculation and 3rd order sensitivities are not implemented.')
            else
                [Gs(iCon), Ds(iCon), Hs(iCon)] = computeObjHess(m_i, con_i, obj_i, opts_i);
            end
        case 'info'
            Fs(iCon) = computeObjInfo(m_i, con_i, obj_i, opts_i, dxdTSol{iCon});
        case 'prob' % Doesn't have a custom computeObj* function
            AbsTol = fixAbsTol(opts_i.AbsTol, 1, opts_i.Continuous, m.nx, 1, opts_i.UseAdjoint, opts_i.UseParams, opts_i.UseSeeds, {opts_i.UseInputControls}, {opts_i.UseDoseControls});
            opts_i.AbsTol = AbsTol{1};
            % Integrate
            ints = integrateAllSys(m, con_i, obj_i, opts_i);
            % Add fields for prior objectives
            [ints.UseParams] = deal(opts_i.UseParams);
            [ints.UseSeeds] = deal(opts_i.UseSeeds);
            [ints.UseInputControls] = deal(opts_i.UseInputControls);
            [ints.UseDoseControls] = deal(opts_i.UseDoseControls);
            % Probability
            p = 1;
            for iObj = 1:length(obj_i)
                p = p * obj_i(iObj).p(ints(iObj));
            end
            ps(iCon) = sum(p);
        case 'logp' % Doesn't have a custom computeObj* function
            AbsTol = fixAbsTol(opts_i.AbsTol, 1, opts_i.Continuous, m.nx, 1, opts_i.UseAdjoint, opts_i.UseParams, opts_i.UseSeeds, {opts_i.UseInputControls}, {opts_i.UseDoseControls});
            opts_i.AbsTol = AbsTol{1};
            % Integrate
            ints = integrateAllSys(m, con_i, obj_i, opts_i);
            % Add fields for prior objectives
            [ints.UseParams] = deal(opts_i.UseParams);
            [ints.UseSeeds] = deal(opts_i.UseSeeds);
            [ints.UseInputControls] = deal(opts_i.UseInputControls);
            [ints.UseDoseControls] = deal(opts_i.UseDoseControls);
            % Probability
            logp = 0;
            for iObj = 1:length(obj_i)
                logp = logp + obj_i(iObj).logp(ints(iObj));
            end
            logps(iCon) = sum(logp);
    end
end

switch mode
    case 'val'
        G = sum(Gs);
        varargout{1} = G;
    case 'grad'
        G = sum(Gs);
        D = opts.fit.Tlocal2T(Ds, 'sum');
        varargout{1} = G;
        varargout{2} = D;
    case 'hess'
        G = sum(Gs);
        D = opts.fit.Tlocal2T(Ds, 'sum');
        H = opts.fit.Tlocal2T(Hs, 'sum');
        varargout{1} = G;
        varargout{2} = D;
        varargout{3} =H;
    case 'info'
        F = opts.fit.Tlocal2T(Fs, 'sum');
        FAll = cell(nCon,1);
        for i = 1:nCon
            FAll{i} = opts.fit.Tlocal2T(Fs(i), 'sum');
        end
        varargout{1} = F;
        varargout{2} = FAll;
    case 'prob'
        p = sum(ps);
        varargout{1} = p;
    case 'logp'
        logp = sum(logps);
        varargout{1} = logp;
end
