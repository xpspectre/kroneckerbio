function [G, D] = computeObjSensAdj(m, con, obj, opts)
verboseAll = max(opts.Verbose-1,0);
if verboseAll; tic; end

%% Set up system
% Constants
nx = m.nx;
ns = m.ns;
nTk = nnz(opts.UseParams);
nTs = nnz(opts.UseSeeds);
nTq = nnz(opts.UseInputControls);
nTh = nnz(opts.UseDoseControls);
nT = nTk + nTs + nTq + nTh;
paramMapper = ParamMapperOneModelType({opts.paramSpec});
nObj = size(obj,1);
d0 = zeros(ns,1); % A fake dose of 0 to subtract out

y = m.y;
dydx = m.dydx;
dydu = m.dydu;
dydk = m.dydk;

% Initialize variables
G = 0;
D = zeros(nT,1);

if opts.Verbose; disp('Integrating adjoint...'); end
if verboseAll; tic; end

% Seeds
s = con.s;

% Input
u = con.u;
d = con.d;

% * Integrate to steady-state
if con.SteadyState
    ssSol = integrateSteadystateSys(m, con, opts);
    ic = ssSol.ye(:,end);
else
    order = 0;
    ic = extractICs(m, con, opts, order);
end

[tF, eve, fin] = collectObservations(m, con, obj(:));

%% Integrate system
% Do not use select methods since the solution is needed at all time
if opts.Continuous
    [der, jac, del] = constructObjectiveSystem();
    sol_sys = accumulateOdeFwdComp(der, jac, 0, tF, [ic; 0], con.Discontinuities, 1:nx, opts.RelTol, opts.AbsTol(1:nx+1), del, eve, fin);
else
    [der, jac, del] = constructSystem();
    sol_sys = accumulateOdeFwdComp(der,jac, 0, tF, ic, con.Discontinuities, 1:nx, opts.RelTol, opts.AbsTol(1:nx), del, eve, fin);
end

% Work down
int_sys = struct;
int_sys.Type = 'Integration.System.Complex';
int_sys.Name = [m.Name ' in ' con.Name];

int_sys.x_names = vec({m.States.Name});
int_sys.u_names = vec({m.Inputs.Name});
int_sys.y_names = vec({m.Outputs.Name});
int_sys.k_names = vec({m.Parameters.Name});
int_sys.s_names = vec({m.Seeds.Name});

int_sys.nx = nx;
int_sys.ny = m.ny;
int_sys.nu = m.nu;
int_sys.nk = m.nk;
int_sys.ns = m.ns;
int_sys.nq = con.nq;
int_sys.nh = con.nh;
int_sys.k = m.k;
int_sys.s = con.s;
int_sys.q = con.q;
int_sys.h = con.h;

int_sys.dydx = m.dydx;
int_sys.dydu = m.dydu;

int_sys.nT = nT;
int_sys.UseParams        = opts.UseParams;
int_sys.UseSeeds         = opts.UseSeeds;
int_sys.UseInputControls = opts.UseInputControls;
int_sys.UseDoseControls  = opts.UseDoseControls;

int_sys.t = sol_sys.x;
int_sys.x = @(t)devals(sol_sys, t);
int_sys.u = con.u;
int_sys.y = @(t)y(t, devals(sol_sys, t), u(t));

int_sys.ie = sol_sys.ie;
int_sys.te = sol_sys.xe;
int_sys.xe = sol_sys.ye;
int_sys.ue = u(int_sys.te);
int_sys.ye = y(int_sys.te, int_sys.xe, int_sys.ue);

int_sys.sol = sol_sys;

% Distribute times for each observation
int_sys = repmat(int_sys, nObj,1);

% Determine which objectives to evaluate for this experiment (those
% that are not objectiveZero)
isobjzero = strcmp('Objective.Data.Zero',{obj(:).Type});
nonZeroObjs = find(~isobjzero);

for j = nonZeroObjs
    if obj(j).Complex
        % Only reveal time points in range of observation
        % Note: deval will still not throw an error outside this range
        int_sys(j).t = [int_sys(j).t(int_sys(j).t < obj(j).tF), obj(j).tF];
    else
        % Evaluate all requested time points
        int_sys(j).t = obj(j).DiscreteTimes;
        int_sys(j).x = int_sys(j).x(int_sys(j).t);
        int_sys(j).u = int_sys(j).u(int_sys(j).t);
        int_sys(j).y = int_sys(j).y(int_sys(j).t);
    end
end

%% Compute G
% Extract continuous term
if opts.Continuous
    G_cont = int_sys(1).sol.y(nx+1,end);
else
    G_cont = 0;
end

% Compute discrete term
G_disc = 0;
discrete_times_all = cell(nObj,1);
for j = nonZeroObjs
    [iDiscG, temp] = obj(j).G(int_sys(j));
    discrete_times_all{j} = row(unique(temp));
    G_disc = G_disc + opts.ObjWeights(j) * iDiscG;
end

discrete_times = vec(unique([discrete_times_all{:}]));

% Add to cumulative goal value
G = G + G_cont + G_disc;

%% Integrate Adjoint
% Construct system
[der, jac, del] = constructAdjointSystem();

% Set initial conditions
ic = zeros(nx+nT,1);

% Integrate [lambda; D] backward in time
sol = accumulateOdeRevSelect(der, jac, 0, tF, ic, [con.Discontinuities; discrete_times], 0, [], opts.RelTol, opts.AbsTol(nx+opts.Continuous+1:nx+opts.Continuous+nx+nT), del);

%% Complete steady-state
if con.SteadyState
    % * Start Adjoint again *
    [der, jac] = constructSteadystateSystem();
    
    % Set initial conditions starting from end of previous run
    ic = sol.y;
    
    % Integrate [lambda; D] backward in time and replace previous run
    sol = accumulateOdeRevSelect(der, jac, 0, ssSol.xe, ic, con.private.BasalDiscontinuities, 0, [], opts.RelTol, opts.AbsTol(nx+opts.Continuous+1:nx+opts.Continuous+nx+nT));
end

%% Add contributions to derivative
% (Subtract contribution because gradient was integrated backward)

% Rate parameters
curD = -sol.y(nx+1:end,end);

% Initial conditions
lambda = -sol.y(1:nx,end);
dx0dk_val = m.dx0dk(s);
curD(1:nTk) = vec(curD(1:nTk)) + dx0dk_val(:,opts.UseParams).' * lambda;
dx0ds_val = m.dx0ds(s);
curD(nTk+1:nTk+nTs) = vec(curD(nTk+1:nTk+nTs)) + dx0ds_val(:,opts.UseSeeds).' * lambda;

% Add to cumulative goal value
D = D + curD;

if opts.Normalized
    T = collectActiveParameters(m, con, opts.UseParams, opts.UseSeeds, {opts.UseInputControls}, {opts.UseDoseControls});
    D = D .* T;
end

if verboseAll; fprintf('|dGdT| = %g\tTime = %0.2f\n', norm(curD), toc); end

if opts.Verbose; fprintf('Summary: |dGdT| = %g\n', norm(D)); end

% Convert gradient to standard form for local T
D = paramMapper.T2Tlocal(D);

% End of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% The system for integrating lambda and D %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [der, jac, del] = constructAdjointSystem()
        dx0dk   = m.dx0dk;
        dx0dd   = m.dx0ds;
        dfdx    = m.dfdx;
        dfdu    = m.dfdu;
        dfdk    = m.dfdk;
        dfdT    = @dfdTSub;
        dudq    = con.dudq;
        dddh    = con.dddh;
        nh      = con.nh;
        
        der = @derivative;
        jac = @jacobian;
        del = @delta;
        
        % Derivative of [lambda; D] with respect to time
        function val = derivative(t, joint)
            u_i = u(t);
            x_i = deval(sol_sys, t, 1:nx);
            y_i = y(t, x_i, u_i);
            dydx_i = dydx(t, x_i, u_i);
            l = joint(1:nx);
            
            % Sum continuous objective functions
            dgdx = zeros(nx,1);
            for i = 1:nObj
                dgdx = dgdx + dydx_i.' * opts.ObjWeights(i)*obj(i).dgdy(t,y_i);
            end
            
            val = [dgdx; zeros(nT,1)] - [dfdx(t,x_i,u_i).'; dfdT(t,x_i,u_i).'] * l;
        end
        
        % Jacobian of [lambda; D] derivative
        function val = jacobian(t, joint)
            ui = u(t);
            x = deval(sol_sys, t, 1:nx);
            
            val = [-dfdx(t,x,ui).', sparse(nx,nT);
                   -dfdT(t,x,ui).', sparse(nT,nT)];
        end
        
        % Discrete effects of the objective function
        function val = delta(t, joint)
            u_i = u(t);
            x_i = deval(sol_sys, t, 1:nx);
            dydx_i = dydx(t, x_i, u_i);
            dydu_i = dydu(t, x_i, u_i);
            dydk_i = dydk(t, x_i, u_i);
            dGdx = zeros(nx,1);
            dGdT = zeros(nT,1);
            for i = 1:nObj
                dGdy = obj(i).dGdy(t, int_sys(i));
                dGdx = dGdx + dydx_i.' * opts.ObjWeights(i)*dGdy;
                dGdk = dydk_i.' * dGdy;
                dGds = zeros(ns, 1);
                dGdq = dudq(t).' * dydu_i.' * dGdy; % q_ % partial dGdq(i)
                dGdh = zeros(nh, 1);
                if t == 0 % Only evaluate the partial derivatives that don't depend on time once at t=0 (t=0 chosen arbitrarily)
                    dGdk = dGdk + obj(i).dGdk(int_sys(i)); % k_ % partial dGdk(i)
                    dGds = dGds + obj(i).dGds(int_sys(i)); % s_ % partial dGds(i)
                    dGdq = dGdq + obj(i).dGdq(int_sys(i)); % q_ % partial dGdq(i)
                    dGdh = dGdh + obj(i).dGdh(int_sys(i)); % h_ % partial dGdh(i)
                end
                dGdT = dGdT + opts.ObjWeights(i)*[dGdk(opts.UseParams); dGds(opts.UseSeeds); dGdq(opts.UseInputControls); dGdh(opts.UseDoseControls)]; % T_ + (k_ -> T_) -> T_
            end
            
            lambda = -joint(1:nx,end) + dGdx; % Update current lambda
            d_i = d(t);
            dx0dk_i = dx0dk(d_i) - dx0dk(d0);
            dx0dk_i = dx0dk_i(:,opts.UseParams);
            dose_change_k = dx0dk_i.' * lambda;
            dx0dd_i = dx0dd(d_i);
            dddh_i = dddh(t);
            dddh_i = dddh_i(:,opts.UseDoseControls);
            dose_change_h = dddh_i.' * dx0dd_i.' * lambda;
            dGdT = dGdT + [dose_change_k; zeros(nTs,1); zeros(nTq,1); dose_change_h];
            
            val = [dGdx; dGdT];
        end
        
        % Modifies dfdk to relate only to the parameters of interest
        function val = dfdTSub(t, x, u)
            val = dfdk(t,x,u);
            dfdq = dfdu(t,x,u) * dudq(t);
            val = [val(:,opts.UseParams), sparse(nx,nTs), dfdq(:,opts.UseInputControls), sparse(nx,nTh)];
        end
    end

    function [der, jac] = constructSteadystateSystem()
        dfdx = m.dfdx;
        dfdu = m.dfdu;
        dfdk = m.dfdk;
        dfdT = @dfdTSub;
        basal_u = con.private.basal_u;
        basal_dudq = con.private.basal_dudq;
        
        der = @derivative;
        jac = @jacobian;
        
        % Derivative of [lambda; D] with respect to time
        function val = derivative(t, joint)
            ui = basal_u(t);
            x = deval(ssSol, t, 1:nx);
            l = joint(1:nx);
            
            val = -[dfdx(-1,x,ui).'; dfdT(t,x,ui).'] * l;
        end
        
        % Jacobian of [lambda; D] derivative
        function val = jacobian(t, joint)
            ui = basal_u(t);
            x = deval(ssSol, t, 1:nx);
            
            val = [-dfdx(-1,x,ui).', sparse(nx,nT);
                   -dfdT(t,x,ui).', sparse(nT,nT)];
        end
        
        % Modifies dfdk to relate only to the parameters of interest
        function val = dfdTSub(t, x, u)
            val = dfdk(-1,x,u);
            dfdq = dfdu(-1,x,u) * basal_dudq(t);
            val = [val(:,opts.UseParams), zeros(nx,nTs), dfdq(:,opts.UseInputControls), sparse(nx,nTh)];
        end
    end

    function [der, jac, del] = constructObjectiveSystem()
        f       = m.f;
        dfdx    = m.dfdx;
        x0      = m.x0;
        
        der = @derivative;
        jac = @jacobian;
        del = @delta;
        
        % Derivative of [x; G] with respect to time
        function val = derivative(t, joint)
            ui = u(t);            
            x = joint(1:nx);
            
            % Sum continuous objective functions
            g = 0;
            for i = 1:n_obj
                g = g + opts.ObjWeights(i) * obj(i).g(t,x,ui);
            end
            
            val = [f(t,x,ui); g];
        end
        
        % Jacobian of [x; G] derivative
        function val = jacobian(t, joint)
            ui = u(t);
            x = joint(1:nx);
            
            % Sum continuous objective gradients
            dgdx = zeros(1,nx);
            for i = 1:n_obj
                dgdx = dgdx + opts.ObjWeights(i) * vec(obj(i).dgdx(t,x,ui)).';
            end
            
            val = [dfdx(t,x,ui), sparse(nx,1);
                          dgdx,            0];
        end

        % Dosing
        function val = delta(t, joint)
            val = [x0(d(t)) - x0(d0); 0];
        end
    end

    function [der, jac, del] = constructSystem()
        f     = m.f;
        dfdx  = m.dfdx;
        x0    = m.x0;
        
        der = @derivative;
        jac = @jacobian;
        del = @delta;
        
        % Derivative of x with respect to time
        function val = derivative(t, x)
            ui   = u(t);
            val = f(t,x,ui);
        end
        
        % Jacobian of x derivative
        function val = jacobian(t, x)
            ui   = u(t);
            val = dfdx(t,x,ui);
        end
        
        % Dosing
        function val = delta(t, x)
            val = x0(d(t)) - x0(d0);
        end
    end
end
