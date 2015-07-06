function [G, D] = computeObjGrad(m, con, obj, opts)
% Compute gradient of objective function w.r.t. params used in condition and
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
AbsTol = fixAbsTol(con.AbsTol, 2, opts.continuous, m.nx, 1, opts.UseAdjoint, opts.UseParams, opts.UseSeeds, {opts.UseInputControls}, {opts.UseDoseControls});
opts.AbsTol = AbsTol{1};
opts.ObjWeights = [obj(:).Weight];

% Switch to Adjoint method if requested
if opts.UseAdjoint
    [G, D] = computeObjSensAdj(m, con, obj, opts);
    return
end

% Continue with forward method

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

if opts.Verbose; fprintf('Integrating sensitivities:\n'); end
verbose_all = max(opts.Verbose-1,0);
if verbose_all; tic; end

% Integrate
ints = integrateAllSens(m, con, obj, opts);

% Extract continuous term
if opts.continuous
    error('Continuous objective functions have not been updated')
    G_cont = ints(1).sol.y(nx+1,end);
    
    dGdTStart = nx+1+nx*nT+1;
    dGdTEnd   = nx+1+nx*nT+nT;
    D_cont = ints(1).sol.y(dGdTStart:dGdTEnd,end);
else
    G_cont = 0;
    D_cont = zeros(nT,1);
end

% Compute discrete term
G_disc = 0;
D_disc = zeros(nT,1);

discrete_times_all = cell(nObj,1);
for i = 1:nObj
    [G_disc_i, temp] = obj(i).G(ints(i));
    discrete_times_all{i} = row(unique(temp));
    G_disc = G_disc + opts.ObjWeights(i) * G_disc_i;
    
    Dk_i = obj(i).dGdk(ints(i));
    Ds_i = obj(i).dGds(ints(i));
    Dq_i = obj(i).dGdq(ints(i));
    Dh_i = obj(i).dGdh(ints(i));
    
    D_disc_i = [Dk_i(opts.UseParams); Ds_i(opts.UseSeeds); Dq_i(opts.UseInputControls); Dh_i(opts.UseDoseControls)];
    
    n_disc = numel(discrete_times_all{i});
    for j = 1:n_disc
        ti = discrete_times_all{i}(j);
        if obj(i).Complex
            dydT_i = reshape(ints(i).dydT(ti), ny, nT);
        else
            dydT_i = reshape(ints(i).dydT(:,j), ny, nT);
        end
        
        dGdy_i = row(obj(i).dGdy(ti, ints(i)));
        
        D_disc_i = D_disc_i + vec(dGdy_i * dydT_i); % _y * y_T -> _T -> T_
    end
    
    D_disc = D_disc + opts.ObjWeights(i) * D_disc_i;
    
end

% Sum discrete and continuous terms
G = G_cont + G_disc;
D = D_cont + D_disc;

if verbose_all; fprintf('|dGdT| = %g\tTime = %0.2f\n', norm(D), toc); end

if opts.Verbose; fprintf('Summary: |dGdT| = %g\n', norm(D)); end

% Convert gradient to standard form for local T
D = paramMapper.T2Tlocal(D);
