function [G, D] = computeObjGrad(m, con, obj, opts)
% Assumes a single m
% Switch to Adjoint method if requested
if opts.UseAdjoint
    warning('UseAdjoint not implemented/modified for changes in fit objective.')
    [G, D] = computeObjSensAdj(m, con, obj, opts);
    return
end

% Continue with forward method
verbose_all = max(opts.Verbose-1,0);

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

% Constants
nx = m.nx;
ny = m.ny;
nTk = nnz(opts.UseParams);
nTs = nnz(opts.UseSeeds);
nTq = nnz(opts.UseInputControls);
nTh = nnz(opts.UseDoseControls);
nT = nTk + nTs + nTq + nTh;
Tmap = {opts.UseParams; opts.UseSeeds; opts.UseInputControls; opts.UseDoseControls};
nObj = length(obj);

if opts.Verbose; fprintf('Integrating sensitivities:\n'); end
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
D = mapT2Tlocal(D, Tmap);

%% Helper function
function Tlocal = mapT2Tlocal(T, Tmap)
% Map combined T vector of fit params for 1 condition to standard form.
% Input:
%   T [ nT x 1 double vector ]
%       Format of combined fit params for this condition
%   Tmap [ 1 x 4 cell vector of nX x 1 logical vectors ]
%       Format of which params are used
% Output:
%   Tlocal [ 1 x 4 cell vector of nX x 1 double vectors ]
%       Format of this condition's fit params
Tlocal = cell(1,4);
Tidx = 0;
for i = 1:4
    Tstart = Tidx + 1;
    Tend = Tstart + nnz(Tmap{i}) - 1;
    Tsub = T(Tstart:Tend);
    
    Tlocal_i = zeros(length(Tmap{i}),1);
    Tlocal_i(Tmap{i}) = Tsub;
    Tlocal{i} = Tlocal_i;
    
    Tidx = Tend;
end