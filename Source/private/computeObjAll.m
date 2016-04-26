function [G, D, H] = computeObjAll(m, con, obj, opts)

nCon = length(con);
ComponentMap = opts.fit.ComponentMap;
[uniqueConInds, ia] = unique(ComponentMap(:,2));
uniqueComponentMap = ComponentMap(ia,:);
assert(nCon == length(uniqueConInds));

if opts.ComputeSensPlusOne
    assert(nargout < 3, 'KroneckerBio:computeObjAll:HessianOrderTooHigh', 'Hessian cannot be computed since ComputeSensPlusOne requires 3rd order sensitivities.')
end

Gs = zeros(nCon,1);
Ds = cell(nCon,1);
Hs = cell(nCon,1);
for iCon = 1:nCon
    m_i = m(uniqueComponentMap(iCon,1));
    con_i = con(iCon);
    obj_i = obj(ComponentMap(:,2) == uniqueComponentMap(iCon,2))';
    
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
    
    if nargout == 1 % objective function value
        if opts.ComputeSensPlusOne
            Gs(iCon) = computeObjGrad(m_i, con_i, obj_i, opts_i);
        else
            Gs(iCon) = computeObj(m_i, con_i, obj_i, opts_i);
        end
    elseif nargout == 2
        if opts.ComputeSensPlusOne
            [Gs(iCon), Ds(iCon)] = computeObjHess(m_i, con_i, obj_i, opts_i);
        else
            [Gs(iCon), Ds(iCon)] = computeObjGrad(m_i, con_i, obj_i, opts_i);
        end
    elseif nargout == 3
        if opts.ComputeSensPlusOne
            error('FitObject:computeObjAll:HessianOrderTooHigh', 'ComputeSensPlusOne = true was specified for an obective Hessian calculation and 3rd order sensitivities are not implemented.')
        else
            [Gs(iCon), Ds(iCon), Hs(iCon)] = computeObjHess(m_i, con_i, obj_i, opts_i);
        end
    end
end

G = sum(Gs);
if nargout > 1
    D = opts.fit.Tlocal2T(Ds, 'sum');
    if nargout > 2
        H = opts.fit.Tlocal2T(Hs, 'sum');
    end
end
