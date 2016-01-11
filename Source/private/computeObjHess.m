function [G, D, H] = computeObjHess(m, con, obj, opts)
% Compute Hessian of objective function w.r.t. params used in condition and
% model
% Inputs:
%   m [ model struct ]
%       Single model
%   con [ condition struct ]
%       Single condition built off model
%   obj [ nObj vector of objective structs ]
%       Objective functions belonging to this condition
%   opts [ struct ]
%       Options struct
%       
% Outputs:
%   G [ scalar double ]
%       Objective function value at params specified in m and con
%   D [ 1 x 1 cell of 4 x 1 cell vector of nXi x 1 double vectors ]
%       Objective function gradient w.r.t. params used in m and con in standard
%       Tlocal form
%   H [ 1 x 1 cell of 4 x 4 cell vector of nXi x nXj double vectors ]
%       Objective function Hessian w.r.t. params used in m and con in standard
%       Tlocal form
%
% Note: this function has been simplified to only allow a single condition

% Process options
% Note: this currently acts as a shim to existing code - a lot of this can be
% cleaned up in the future
opts.UseParams        = logical(con.ParamsSpec{1});
opts.UseSeeds         = logical(con.ParamsSpec{2});
opts.UseInputControls = logical(con.ParamsSpec{3});
opts.UseDoseControls  = logical(con.ParamsSpec{4});

opts.continuous = any(obj(:).Continuous);
opts.RelTol = con.RelTol;
AbsTol = fixAbsTol(con.AbsTol, 3, opts.continuous, m.nx, 1, opts.UseAdjoint, opts.UseParams, opts.UseSeeds, {opts.UseInputControls}, {opts.UseDoseControls});
opts.AbsTol = AbsTol{1};
opts.ObjWeights = [obj(:).Weight];

% Constants
nx = m.nx;
ny = m.ny;
nTk = nnz(opts.UseParams);
nTs = nnz(opts.UseSeeds);
nTq = nnz(opts.UseInputControls);
nTh = nnz(opts.UseDoseControls);
nT = nTk + nTs + nTq + nTh;
paramMapper = ParamMapperOneModelType(con);
nObj = length(obj);

if opts.Verbose; fprintf('Integrating curvature:\n'); end
verbose_all = max(opts.Verbose-1,0);
if verbose_all; tic; end

% Integrate
ints = integrateAllCurv(m, con, obj, opts);

% NLME-specific fields
if is(m, 'Model.Nlme')
    ints.fOmega = m.fOmega;
end

% Extract continuous term
if opts.continuous
    error('Continuous objective functions have not been updated')
    G_cont = ints(1).sol.y(nx+1,end);
    
    dGdTStart = nx+1+nx*nT+1;
    dGdTEnd   = nx+1+nx*nT+nT;
    D_cont = ints(1).sol.y(dGdTStart:dGdTEnd,end);
    
    d2GdT2Start = nx+1+nx*nT+nx*nT*nT+1;
    d2GdT2End   = nx+1+nx*nT+nx*nT*nT+nT;
    H_cont = reshape(ints(1).sol.y(d2GdT2Start:d2GdT2End,end), nT,nT);
else
    G_cont = 0;
    D_cont = zeros(nT,1);
    H_cont = zeros(nT,nT);
end

% Compute discrete term
G_disc = 0;
D_disc = zeros(nT,1); % T_
H_disc = zeros(nT,nT); % T_T

discrete_times_all = cell(nObj,1);
for i = 1:nObj
    [G_disc_i, temp] = obj(i).G(ints(i));
    discrete_times_all{i} = row(unique(temp));
    G_disc = G_disc + opts.ObjWeights(i) * G_disc_i;
    
    Dk_i = obj(i).dGdk(ints(i));
    Ds_i = obj(i).dGds(ints(i));
    Dq_i = obj(i).dGdq(ints(i));
    Dh_i = obj(i).dGdh(ints(i));
    if opts.Normalized
        Dk_i = Dk_i.*ints(i).k;
        Ds_i = Ds_i.*ints(i).s;
        Dq_i = Dq_i.*ints(i).q;
        Dh_i = Dh_i.*ints(i).h;
    end
    D_disc_i = [Dk_i(opts.UseParams); Ds_i(opts.UseSeeds); Dq_i(opts.UseInputControls); Dh_i(opts.UseDoseControls)];
    
    Hk_i = obj(i).d2Gdk2(ints(i));
    Hs_i = obj(i).d2Gds2(ints(i));
    Hq_i = obj(i).d2Gdq2(ints(i));
    Hh_i = obj(i).d2Gdh2(ints(i));
    if opts.Normalized
        Hk_i = spdiags(ints(i).k,0,m.nk,m.nk) * Hk_i * spdiags(ints(i).k,0,m.nk,m.nk) + spdiags(Dk_i,0,m.nk,m.nk);
        Hs_i = spdiags(ints(i).s,0,m.ns,m.ns) * Hs_i * spdiags(ints(i).s,0,m.ns,m.ns) + spdiags(Ds_i,0,m.ns,m.ns);
        Hq_i = spdiags(ints(i).q,0,con.nq,con.nq) * Hq_i * spdiags(ints(i).q,0,con.nq,con.nq) + spdiags(Dq_i,0,con.nq,con.nq);
        Hh_i = spdiags(ints(i).h,0,con.nh,con.nh) * Hh_i * spdiags(ints(i).h,0,con.nh,con.nh) + spdiags(Dh_i,0,con.nh,con.nh);
    end
    H_disc_i = [Hk_i(opts.UseParams,opts.UseParams), zeros(nTk,nT-nTk);
        zeros(nTs,nTk), Hs_i(opts.UseSeeds,opts.UseSeeds), zeros(nTs,nTq+nTh);
        zeros(nTq,nTk+nTs), Hq_i(opts.UseInputControls,opts.UseInputControls), zeros(nTq,nTh);
        zeros(nTh,nT-nTh), Hh_i(opts.UseDoseControls,opts.UseDoseControls)];
    
    n_disc = numel(discrete_times_all{i});
    for j = 1:n_disc
        ti = discrete_times_all{i}(j);
        if obj(i).Complex
            dydT_i = reshape(ints(i).dydT(ti), ny, nT);
            d2ydT2_i = reshape(ints(i).d2ydT2(ti), ny, nT*nT);
        else
            dydT_i = reshape(ints(i).dydT(:,j), ny, nT);
            d2ydT2_i = reshape(ints(i).d2ydT2(:,j), ny, nT*nT);
        end
        
        dGdy_i = row(obj(i).dGdy(ti, ints(i)));
        d2Gdy2_i = obj(i).d2Gdy2(ti, ints(i));
        
        % D = sum(dGdy(t) *{y.y} dydT(t), t) + dGdT
        D_disc_i = D_disc_i + vec(dGdy_i * dydT_i); % _y * y_T -> _T -> T_
        
        % H = (d2Gdy2(t) *{y.y} dydT(t) *{y.y}) dydT(t) + dGdy *{y.y} d2ydT2 + d2GdT2
        H_disc_i = H_disc_i + dydT_i.' * d2Gdy2_i * dydT_i + reshape(dGdy_i * d2ydT2_i, nT, nT); % (y_T -> T_y) * y_y * y_T -> T_T
    end
    
    D_disc = D_disc + opts.ObjWeights(i) * D_disc_i;
    H_disc = H_disc + opts.ObjWeights(i) * H_disc_i;
end

% Sum discrete and continuous terms
G = G_cont + G_disc;
D = D_cont + D_disc;
H = H_cont + H_disc;

if verbose_all; fprintf('|dGdT| = %g\t||d2GdT2|| = %g\tTime = %0.2f\n', norm(D), det(H), toc); end

if opts.Verbose; fprintf('Summary: |dGdT| = %g\t||d2GdT2|| = %g\n', norm(D), det(H)); end

% Convert gradient and Hessian to standard form for local T
D = paramMapper.T2Tlocal(D);
H = paramMapper.T2Tlocal(H);