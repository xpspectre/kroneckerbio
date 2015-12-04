function obs = observationFOCEI(output, timelist, method, name)
% Nonlinear mixed effects observation and objective function for a patient
%
% Inputs:
%   output [ string | cell array of strings ]
%       Name of output{s} of interest
%   timelist [ nt x 1 vector of nonnegative doubles ]
%       Vector of times of interest
%   method [ {'FO'} | 'FOCE' ]
%       Approximation method
%   name [ string {'NLME'} ]
%       String name of model
% Outputs:
%   obj [ Objective.NLME struct ]
%       Observation struct
%
% Note/TODO: Currently only supports a single output. Extend to use an arbitrary
% number of outputs later.

% Clean up inputs
if nargin < 4
    name = [];
    if nargin < 3
        method = [];
    end
end

if isempty(method)
    method = 'FOCEI';
end
if isempty(name)
    name = 'FOCEI';
end

% Clean up outputs
output = fixOutputName(output);
nOutputs = length(output);

% Find unique timelist
discrete_times = row(unique(timelist));

obs = [];
obs.Type    = 'Observation.Data.FOCEI';
obs.Name    = name;
obs.Complex = false;

obs.tF            = max([0, discrete_times]);
obs.DiscreteTimes = discrete_times;
obs.Outputs       = output;

obs.Simulation = @simulation;
obs.Objective  = @objective;

obs = pastestruct(observationZero(), obs);

    function sim = simulation(int)
        sim.Type = 'Simulation.Data.FOCEI';
        sim.Name = name;
        sim.Output = output;
        sim.Times = discrete_times;
        
        y_Ind = lookupOutput(output, int); % take actual output being compared to measurements (w/ error model from eps)
        sim.Predictions = int.y(y_Ind,:);
    end

    function obj = objective(measurements)
        obj = objectiveNLME(output, timelist, measurements, method, name);
    end

end

function obj = objectiveNLME(output, timelist, measurements, method, name)
% Inputs:
%   output
%   timelist
%   measurements [ ni x nOutpus matrix of doubles ]

% Clean up outputs
output = fixOutputName(output);
nOutputs = length(output);

% Find unique timelist
discrete_times = row(unique(timelist));

% Verify measurements match times
assert(all(size(measurements) == [length(timelist), nOutputs]), 'objectiveFOCEI:invalidMeasurements', 'Dimensions of measurements don''t match nTimes and/or nOutputs')

% Inherit observation
obj = observationFOCEI(output, timelist, method, name);

obj.Type = 'Objective.Data.NLME';
obj.Continuous = false;

obj.G      = @G;
obj.dGdk   = @dGdk;
obj.d2Gdk2 = @d2Gdk2;

obj.p    = @p;
obj.logp = @logp;
obj.F    = @F;
obj.Fn   = @Fn;

obj = pastestruct(objectiveZero(), obj);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Parameter fitting functions %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Debug: offload all calcs, of different orders, to this function
% Notes:
%   only use squeeze when final result is 1 or 2 dimensional
%       using reshape can remove singleton dimensions unnecessarily when
%       something has 1 dim but could have more in other cases, like scalar
%       versus matrix Rij
%   partial derivatives from function handles; total derivative from filled in
%   values
%       partial derivatives start with pd...
%       filling in values for partial derivatives denoted by starting with p
%   assume Rij is diagonal
%       but maintain it as a matrix for various operations
%   dOmega_dtheta = 0 assumed for now
    function [val, discrete] = calcParts(int, order)
        % All terms with star are evaluated at optimal eta star (same as eta in FO method)
        % Inputs:
        %   int [ struct ]
        %       Integrated solution struct
        %   order [ 1 | 2 ]
        %       Whether to require sensitivity (1) or curvature (2). Redundant
        %       check with presence of int.d2ydx2
        
        % Debug: display current order of int
        if isfield(int, 'd2ydx2')
            fprintf('Calculating obj fun and grad with 2nd order sensitivities\n')
        else
            fprintf('Calculating obj fun and grad with 1st order sensitivities\n')
        end
        
        %         m = int.m; % may just add the handles to int directly instead of models
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get common components
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        nx = int.nx;
        ny = int.ny;
        nT = int.nT;
        
        % Indexing
        [RiInds, hInds, errInds] = getOutputInds(int.y_names);
        nh = length(hInds);
        [thetaInds, etaInds, omegaInds, sigmaInds] = getNLMEInds(int.k_names); % no more eps
        ntheta = length(thetaInds);
        neta = length(etaInds);
        
        Omega = getOmega(int.k, int.k_names, thetaInds, omegaInds);
        Omega = Omega.^2; % TODO: variances are omegas squared?
        Omega_ = inv(Omega);
        
        etai = int.k(etaInds);
        
        % Do everything by timepoints
        ni = size(measurements,1);
        
        if order >= 1
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Assemble G components
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            epsi        = zeros(nh,ni);
            
            Ri          = zeros(nh,nh,ni);
            Ri_         = Ri;
            
            pdh_detai   = zeros(nh,neta,ni);
            pdh_dxi     = zeros(nh,nx,ni);
            pdRi_detai  = zeros(nh,nh,neta,ni);
            pdRi_dxi    = zeros(nh,nh,nx,ni);
            dxi_detai   = zeros(nx,neta,ni);
            
            depsi_detai = zeros(nh,neta,ni);
            dRi_detai   = zeros(nh,nh,neta,ni);
            
            a           = zeros(nh,neta,ni);
            B           = zeros(nh,nh,ni);
            c           = zeros(nh,nh,neta,ni);
            
            for j = 1:ni %%%DEBUG%%%
                
                % Eq 5
                %   Measurements is a ni x nh matrix (transposed from everything else)
                %   yhat used here is E[y], excludes intra-individual error term
                epsij = zeros(nh,1);
                for ih = 1:nh
                    epsij(ih) = measurements(j,ih) - int.y(hInds(ih),j);
                end
                epsi(:,j) = epsij;
                
                % Partial derivatives from function handles, time dependence from
                %   specifying x and u at time indices
                pdydxj = int.dydx([], int.x(:,j), int.u(:,j));
                pdydkj = int.dydk([], int.x(:,j), int.u(:,j));
                
                % Get Ri, intra-individual variability model covariance matrix
                %   Assume Rij is diagonal
                Rij = zeros(nh);
                for ih = 1:nh
                    Rij(ih,ih) = int.y(RiInds(ih),j);
                end
                Ri(:,:,j) = Rij;
                
                % Get inverse of Ri once for convenience
                Rij_ = inv(Rij);
                Ri_(:,:,j) = Rij_;
                
                % Eq 16
                B(:,:,j) = 2*Rij_;
                
                % In Eq 19
                pdh_detaik = zeros(nh,neta);
                for ih = 1:nh
                    for k = 1:neta
                        pdh_detaik(ih,k) = pdydkj(hInds(ih),etaInds(k));
                    end
                end
                pdh_detai(:,:,j) = pdh_detaik;
                
                % In Eq 19
                pdh_dxij = zeros(nh,nx);
                for ih = 1:nh
                    for ix = 1:nx
                        pdh_dxij(ih,ix) = pdydxj(hInds(ih),ix);
                    end
                end
                pdh_dxi(:,:,j) = pdh_dxij;
                
                % In Eq 20
                pdRij_detaik = zeros(nh,nh,neta);
                for ih = 1:nh
                    for k = 1:neta
                        pdRij_detaik(ih,ih,k) = pdydkj(RiInds(ih),etaInds(k));
                    end
                end
                pdRi_detai(:,:,:,j) = pdRij_detaik;
                
                % In Eq 20
                pdRij_dxij = zeros(nh,nh,nx);
                for ih = 1:nh
                    for ix = 1:nx
                        pdRij_dxij(ih,ih,ix) = pdydxj(RiInds(ih),ix);
                    end
                end
                pdRi_dxi(:,:,:,j) = pdRij_dxij;
                
                % In Eq 19,20
                dxij_detaik = zeros(nx,neta);
                for ix = 1:nx
                    for k = 1:neta
                        ind = nx*(etaInds(k)-1) + ix;
                        dxij_detaik(ix,k) = int.dxdT(ind,j);
                    end
                end
                dxi_detai(:,:,j) = dxij_detaik;
                
            end
            
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Calculate G
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for j = 1:ni %%%DEBUG%%%
                for k = 1:neta
                    
                    % Eq 19
                    pdh_detaik  = squeeze(pdh_detai(:,k,j));
                    pdh_dxij    = squeeze(pdh_dxi(:,:,j));
                    dxij_detaik = squeeze(dxi_detai(:,k,j));
                    depsij_detaik = -(pdh_detaik + pdh_dxij*dxij_detaik);
                    depsi_detai(:,k,j) = depsij_detaik;
                    
                    % Eq 20
                    % Careful with squeezing/reshaping Ri with single output when tensor is needed - will remove extra singleton dims
                    pdRij_detaik = squeeze(pdRi_detai(:,:,k,j));
                    pdRij_dxij   = reshape(pdRi_dxi(:,:,:,j), nh,nh,nx);
                    dRij_detaik = pdRij_detaik + tprod(pdRij_dxij, [1,2,-1], dxij_detaik, [-1]); % contract on xi
                    dRi_detai(:,:,k,j) = dRij_detaik;
                    
                    % Eq 15
                    epsij = squeeze(epsi(:,j));
                    Rij_  = squeeze(Ri_(:,:,j));
                    a(:,k,j) = depsij_detaik' - epsij'*Rij_*dRij_detaik;
                    
                    % Eq 17
                    c(:,:,k,j) = Rij_*dRij_detaik;
                    
                end
            end
            
            % Eq 14, Almquist Hessian approx
            Hi = zeros(neta);
            for k = 1:neta
                for l = 1:neta
                    Hikl = zeros(ni,1);
                    for j = 1:ni
                        al = squeeze(a(:,l,j));
                        ak = squeeze(a(:,k,j));
                        Bj = squeeze(B(:,:,j));
                        cl = squeeze(c(:,:,l,j));
                        ck = squeeze(c(:,:,k,j));
                        Hikl(j) = al*Bj*ak' + trace(-cl*ck);
                    end
                    Hi(k,l) = -1/2 * sum(Hikl) - Omega_(k,l);
                end
            end
            
            % Eq 70, NONMEM Hessian approx (simpler)
            Hi_tilde = zeros(neta);
            for k = 1:neta
                for l = 1:neta
                    Hikl = zeros(ni,1);
                    for j = 1:ni
                        depsij_detail = squeeze(depsi_detai(:,l,j));
                        depsij_detaik = squeeze(depsi_detai(:,k,j));
                        cl            = squeeze(c(:,:,l,j));
                        ck            = squeeze(c(:,:,k,j));
                        Rij_          = squeeze(Ri_(:,:,j));
                        Hikl(j) = depsij_detail'*Rij_*depsij_detaik + 1/2 * trace(cl*ck);
                    end
                    Hi_tilde(k,l) = -sum(Hikl) - Omega_(k,l);
                end
            end
            
            % Eq 10
            lij = zeros(ni,1);
            for j = 1:ni
                epsij = squeeze(epsi(:,j));
                Rij   = squeeze(Ri(:,:,j));
                Rij_  = squeeze(Ri_(:,:,j));
                lij(j) = epsij'*Rij_*epsij + log(det(2*pi*Rij));
            end
            li = -1/2 * sum(lij) - 1/2 * etai'*Omega_*etai - 1/2 * log(det(2*pi*Omega));
            
            % Eq 18, for each patient
            logLFi = li - 1/2 * log(det(-Hi/(2*pi))); % Almquist
            %             logLFi = li - 1/2 * log(det(-Hi_tilde/(2*pi))); % NONMEM
            
            val = logLFi
            %             val = -2*logLFi % NONMEM minimizes -2loglikelihood
            discrete = 0;
        end
        
        % Reevaluate int at etastar vs eta in FOCE vs FO
        
        if order >= 2
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Assemble dGdk components
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Assume
            Ristar  = Ri;
            Ristar_ = Ri_;
            Rijstar = Rij;
            Rijstar_ = Rij_;
            dOmega_dtheta = zeros(neta,neta,ntheta); % assume
            dOmega_dtheta_ = zeros(neta,neta,ntheta); % assume (calculated from inv of Omega, then set to 0)
            
            pdh_dtheta = zeros(nh,ntheta,ni); % depends on time
            pdRi_dtheta = zeros(nh,nh,ntheta,ni);
            dxi_dtheta = zeros(nx,ntheta,ni);
            
            pd2h_detaidtheta = zeros(nh,neta,ntheta,ni);
            pd2Ri_detaidtheta = zeros(nh,nh,neta,ntheta,ni);
            
            pd2h_detaidxi = zeros(nh,neta,nx,ni);
            pd2Ri_detaidxi = zeros(nh,nh,neta,nx,ni);
            
            pd2h_dxidtheta = zeros(nh,nx,ntheta,ni);
            pd2Ri_dxidtheta = zeros(nh,nh,nx,ntheta,ni);
            
            pd2h_dxi2 = zeros(nh,nx,nx,ni);
            pd2Ri_dxi2 = zeros(nh,nh,nx,nx,ni);
            
            pd2h_detaidetai = zeros(nh,neta,neta,ni); % depends on time
            pd2Ri_detaidetai = zeros(nh,nh,neta,neta,ni);
            
            pd2h_dxidetai = zeros(nh,nx,neta,ni);
            pd2Ri_dxidetai = zeros(nh,nh,nx,neta,ni);
            
            d2xi_detaidtheta = zeros(nx,neta,ntheta,ni);
            d2xi_detaidetai = zeros(nx,neta,neta,ni);
            
            for j = 1:ni
                
                
                % In Eq 34
                pdh_dthetaj = zeros(nh,ntheta);
                for ih = 1:nh
                    for m = 1:ntheta
                        pdh_dthetaj(ih,m) = pdydkj(hInds(ih),thetaInds(m));
                    end
                end
                pdh_dtheta(:,:,j) = pdh_dthetaj;
                
                % In Eq 35
                pdRij_dtheta = zeros(nh,nh,ntheta);
                for ih = 1:nh
                    for m = 1:ntheta
                        pdRij_dtheta(ih,ih,m) = pdydkj(hInds(ih),thetaInds(m));
                    end
                end
                pdRi_dtheta(:,:,:,j) = pdRij_dtheta;
                
                % In Eq 34
                dxij_dtheta = zeros(nx,ntheta);
                for ix = 1:nx
                    for m = 1:ntheta
                        ind = nx*(thetaInds(m)-1) + ix;
                        dxij_dtheta(ix,m) = int.dxdT(ind,j);
                    end
                end
                dxi_dtheta(:,:,j) = dxij_dtheta;
                
                % 2nd partial derivatives from function handles, time dependence from
                %   specifying x and u at time indices
                % Ordering, as expected from order of derivatives (denom right to left):
                %   1st dim: num and 2nd denom(left), num changes faster
                %   2nd dim: 1st denom(right)
                % In the extracted "Hessian":
                %   1st dim: left denom
                %   2nd dim: right denom
                pd2ydx2j  = int.d2ydx2([], int.x(:,j), int.u(:,j)); % [[y,x2], x1]
                pd2ydk2j  = int.d2ydk2([], int.x(:,j), int.u(:,j)); % [[y,k2], k1]
                pd2ydkdxj = int.d2ydkdx([], int.x(:,j), int.u(:,j)); % [[y,x], k]
                pd2ydxdkj = int.d2ydxdk([], int.x(:,j), int.u(:,j)); % [[y,k], x]
                
                % In Eq 37
                pd2h_detaidtheta   = zeros(nh,neta,ntheta);
                pd2Rij_detaidtheta = zeros(nh,nh,neta,ntheta);
                for ih = 1:nh
                    for k = 1:neta
                        for m = 1:ntheta
                            etaInd = etaInds(k);
                            hthetaInd = ny*(thetaInds(m)-1) + hInds(ih);
                            RithetaInd = ny*(thetaInds(m)-1) + RiInds(ih);
                            
                            pd2h_detaidtheta(ih,k,m) = pd2ydk2j(hthetaInd,etaInd);
                            pd2Rij_detaidtheta(ih,ih,k,m) = pd2ydk2j(RithetaInd,etaInd);
                        end
                    end
                end
                pd2h_detaidtheta(:,:,:,j) = pd2h_detaidtheta;
                pd2Ri_detaidtheta(:,:,:,:,j) = pd2Rij_detaidtheta;
                
                % In Eq 37
                pd2h_detaidxij = zeros(nh,neta,nx);
                pd2Rij_detaidxij = zeros(nh,nh,neta,nx);
                for ih = 1:nh
                    for k = 1:neta
                        for ix = 1:nx
                            etaInd = etaInds(k);
                            hxInd = ny*(ix-1) + hInds(ih);
                            RixInd = ny*(ix-1) + RiInds(ih);
                            
                            pd2h_detaidxij(ih,k,ix) = pd2ydkdxj(hxInd,etaInd);
                            pd2Rij_detaidxij(ih,ih,k,ix) = pd2ydkdxj(RixInd,etaInd);
                        end
                    end
                end
                pd2h_detaidxi(:,:,:,j) = pd2h_detaidxij;
                pd2Ri_detaidxi(:,:,:,:,j) = pd2Rij_detaidxij;
                
                % In Eq 37
                pd2h_dxijdtheta = zeros(nh,nx,ntheta);
                pd2Rij_dxijdtheta = zeros(nh,nh,nx,ntheta);
                for ih = 1:nh
                    for ix = 1:nx
                        for m = 1:ntheta
                            xInd = ix;
                            hthetaInd = ny*(thetaInds(m)-1) + hInds(ih);
                            RithetaInd = ny*(thetaInds(m)-1) + RiInds(ih);
                            
                            pd2h_dxijdtheta(nh,nx,ntheta) = pd2ydxdkj(hthetaInd,xInd);
                            pd2Rij_dxijdtheta(nh,nh,nx,ntheta) = pd2ydxdkj(RithetaInd,xInd);
                        end
                    end
                end
                pd2h_dxidtheta(:,:,:,j) = pd2h_dxijdtheta;
                pd2Ri_dxidtheta(:,:,:,:,j) = pd2Rij_dxijdtheta;
                
                % In Eq 37
                pd2h_dxij2 = zeros(nh,nx,nx);
                pd2Rij_dxij2 = zeros(nh,nh,nx,nx);
                for ih = 1:nh
                    for ix = 1:nx
                        for jx = 1:nx
                            xInd = ix;
                            hxInd = ny*(jx-1) + hInds(ih);
                            RixInd = ny*(jx-1) + RiInds(ih);
                            
                            pd2h_dxij2(ih,ix,jx) = pd2ydx2j(hxInd,xInd);
                            pd2Rij_dxij2(ih,ih,ix,jx) = pd2ydx2j(RixInd,xInd);
                        end
                    end
                end
                pd2h_dxi2(:,:,:,j) = pd2h_dxij2;
                pd2Ri_dxi2(:,:,:,:,j) = pd2Rij_dxij2;
                
                % In Eq 37
                pd2h_detaidetaij = zeros(nh,neta,neta); % depends on time
                pd2Rij_detaidetai = zeros(nh,nh,neta,neta);
                for ih = 1:nh
                    for k = 1:neta
                        for l = 1:neta
                            etaInd = etaInds(k);
                            hetaInd = ny*(etaInds(l)-1) + hInds(ih);
                            RietaInd = ny*(etaInds(l)-1) + RiInds(ih);
                            
                            pd2h_detaidetaij(ih,k,l) = pd2ydk2j(hetaInd,etaInd);
                            pd2Rij_detaidetai(ih,ih,k,l) = pd2ydk2j(RietaInd,etaInd);
                        end
                    end
                end
                pd2h_detaidetai(:,:,:,j) = pd2h_detaidetaij;
                pd2Ri_detaidetai(:,:,:,:,j) = pd2Rij_detaidetai;
                
                % In Eq 37
                pd2h_dxijdetai = zeros(nh,nx,neta);
                pd2Rij_dxijdetai = zeros(nh,nh,nx,neta);
                for ih = 1:nh
                    for ix = 1:nx
                        for l = 1:neta
                            xInd = ix;
                            hetaInd = ny*(etaInds(l)-1) + hInds(ih);
                            RietaInd = ny*(etaInds(l)-1) + RiInds(ih);
                            
                            pd2h_dxijdetai(ih,ix,l) = pd2ydxdkj(hetaInd,xInd);
                            pd2Rij_dxijdetai(ih,ih,ix,l) = pd2ydxdkj(RietaInd,xInd);
                        end
                    end
                end
                pd2h_dxidetai(:,:,:,j) = pd2h_dxijdetai;
                pd2Ri_dxidetai(:,:,:,:,j) = pd2Rij_dxijdetai;
                
                % In Eq 37
                % x changes fastest, then thetam; etaik changes slowest
                d2xij_detaikdthetam = zeros(nx,neta,ntheta);
                for ix = 1:nx
                    for k = 1:neta
                        for m = 1:ntheta
                            ind = nx*nT*(etaInds(k)-1) + nx*(thetaInds(m)-1) + ix;
                            d2xij_detaikdthetam(ix,k,m) = int.d2xdT2(ind,j); % 2 x 13 x 13 x time; 1st 3 dims compressed into row indices
                        end
                    end
                end
                d2xi_detaidtheta(:,:,:,j) = d2xij_detaikdthetam;
                
                % In Eq 37
                % x changes fastest, then etail; etaik changes slowest
                d2xij_detaikdetail = zeros(nx,neta,neta);
                for ix = 1:nx
                    for k = 1:neta
                        for l = 1:neta
                            ind = nx*nT*(etaInds(k)-1) + nx*(etaInds(l)-1) + ix;
                            d2xij_detaikdetail(ix,k,l) = int.d2xdT2(ind,j); % 2 x 13 x 13 x time; 1st 3 dims compressed into row indices
                        end
                    end
                end
                d2xi_detaidetai(:,:,:,j) = d2xij_detaikdetail;
            end
            
            
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Calculate dGdk
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Assume
            epsistar = epsi;
            astar = a;
            Bstar = B;
            cstar = c;
            
            depsi_dtheta = zeros(nh,ntheta,ni);
            dRi_dtheta   = zeros(nh,nh,ntheta,ni);
            for j = 1:ni
                for m = 1:ntheta
                    
                    % Eq 34
                    pdh_dthetam  = squeeze(pdh_dtheta(:,m));
                    pdh_dxij     = squeeze(pdh_dxi(:,:,j));
                    dxij_dthetam = squeeze(dxi_dtheta(:,m,j));
                    depsij_dthetam = -(pdh_dthetam + pdh_dxij*dxij_dthetam);
                    depsi_dtheta(:,m,j) = depsij_dthetam;
                    
                    % Eq 36
                    pdRij_dthetam = squeeze(pdRi_dtheta(:,:,m,j));
                    pdRij_dxij    = reshape(pdRi_dxi(:,:,:,j), nh,nh,nx);
                    dRij_dthetam = pdRij_dthetam + tprod(pdRij_dxij, [1,2,-1], dxij_dthetam, [-1]);
                    dRi_dtheta(:,:,m,j) = dRij_dthetam;
                end
            end
            
            % Eq 46 second derivs for Eq 13
            % Eq 19 deriv wrt eta again
            % Eq 20 deriv wrt eta again
            d2epsi_detaidetai = zeros(nh,neta,neta,ni);
            d2Ri_detaidetai   = zeros(nh,nh,neta,neta,ni);
            for j = 1:ni
                
                for l = 1:neta
                    for k = 1:neta
                        dxij_detaik        = squeeze(dxi_detai(:,k,j));
                        dxij_detail        = squeeze(dxi_detai(:,l,j));
                        d2xij_detaikdetail = squeeze(d2xi_detaidetai(:,k,l,j));
                        
                        pdh_dxij          = squeeze(pdh_dxi(:,:,j));
                        pd2h_dxij2        = reshape(pd2h_dxi2(:,:,:,j), nh,nx,nx);
                        pd2h_dxijdetail   = squeeze(pd2h_dxidetai(:,:,l,j));
                        pd2h_detaikdxij   = reshape(pd2h_detaidxi(:,k,:,j), nh,nx);
                        pd2h_detaikdetail = squeeze(pd2h_detaidetai(:,k,l,j));
                        
                        d2epsi_detaidetai(:,k,l,j) = -(pd2h_detaikdetail + pd2h_detaikdxij*dxij_detail + ...
                            (pd2h_dxijdetail + tprod(pd2h_dxij2, [1,2,-1], dxij_detail, [-1]))*dxij_detaik + ...
                            pdh_dxij*d2xij_detaikdetail);
                        
                        pdRij_dxij          = reshape(pdRi_dxi(:,:,:,j), nh,nh,nx);
                        pd2Rij_dxij2        = reshape(pd2Ri_dxi2(:,:,:,:,j), nh,nh,nx,nx);
                        pd2Rij_dxijdetail   = reshape(pd2Ri_dxidetai(:,:,:,l,j), nh,nh,nx);
                        pd2Rij_detaikdxij   = reshape(pd2Ri_detaidxi(:,:,k,:,j), nh,nh,nx);
                        pd2Rij_detaikdetail = squeeze(pd2Ri_detaidetai(:,:,k,l,j));
                        
                        d2Ri_detaidetai(:,:,k,l,j) = pd2Rij_detaikdetail + tprod(pd2Rij_detaikdxij, [1,2,-1], dxij_detail, [-1]) + ...
                            tprod((pd2Rij_dxijdetail + tprod(pd2Rij_dxij2, [1,2,3,-1], dxij_detail, [-1])), [1,2,-1], dxij_detaik, [-1]) + ...
                            tprod(pdRij_dxij, [1,2,-1], d2xij_detaikdetail, [-1]);
                    end
                end
            end
            
            % For Eq 37 term1 and Eq 47
            % For Eq 37 corresponding Rij term1 and Eq 47
            d2epsi_detaidtheta = zeros(nh,neta,ntheta,ni);
            d2Ri_detaidtheta = zeros(nh,nh,neta,ntheta,ni);
            for j = 1:ni
                for m = 1:ntheta
                    for k = 1:neta
                        dxij_detaik         = squeeze(dxi_detai(:,k,j));
                        dxij_dthetam        = squeeze(dxi_dtheta(:,m,j));
                        d2xij_detaikdthetam = squeeze(d2xi_detaidtheta(:,k,m,j));
                        
                        pdh_dxij            = squeeze(pdh_dxi(:,:,j));
                        pd2h_dxij2          = reshape(pd2h_dxi2(:,:,:,j), nh,nx,nx);
                        pd2h_dxijdthetam    = squeeze(pd2h_dxidtheta(:,:,m,j));
                        pd2h_detaikdxij     = reshape(pd2h_detaidxi(:,k,:,j), nh,nx);
                        pd2h_detaikdthetam  = squeeze(pd2h_detaidtheta(:,k,m,j));
                        
                        d2epsi_detaidtheta(:,k,m,j) = -(pd2h_detaikdthetam + pd2h_detaikdxij*dxij_dthetam + ...
                            (pd2h_dxijdthetam + tprod(pd2h_dxij2, [1,2,-1], dxij_dthetam, [-1]))*dxij_detaik + ...
                            pdh_dxij*d2xij_detaikdthetam);
                        
                        pdRij_dxij           = reshape(pdRi_dxi(:,:,:,j), nh,nh,nx);
                        pd2Rij_dxij2         = reshape(pd2Ri_dxi2(:,:,:,:,j), nh,nh,nx,nx);
                        pd2Rij_dxijdthetam   = squeeze(pd2Ri_dxidtheta(:,:,nx,m,j));
                        pd2Rij_detaikdxij    = reshape(pd2Ri_detaidxi(:,:,k,:,j), nh,nh,nx);
                        pd2Rij_detaikdthetam = squeeze(pd2Ri_detaidtheta(:,:,k,m,j));
                        
                        d2Ri_detaidtheta(:,:,k,m,j) = pd2Rij_detaikdthetam + tprod(pd2Rij_detaikdxij, [1,2,-1], dxij_dthetam, [-1]) + ...
                            tprod((pd2Rij_dxijdthetam + tprod(pd2Rij_dxij2, [1,2,3,-1], dxij_dthetam, [-1])), [1,2,-1], dxij_detaik, [-1]) + ...
                            tprod(pdRij_dxij, [1,2,-1], d2xij_detaikdthetam, [-1]);
                    end
                end
            end
            
            % Eq 13
            d2li_detai2 = zeros(neta,neta);
            for k = 1:neta
                for l = 1:neta
                    
                    terms = zeros(ni,1);
                    for j = 1:ni
                        epsij                = squeeze(epsi(:,j));
                        Rij_                 = squeeze(Ri_(:,:,j));
                        depsij_detail        = squeeze(depsi_detai(:,l,j));
                        depsij_detaik        = squeeze(depsi_detai(:,k,j));
                        dRij_detail          = squeeze(dRi_detai(:,:,l,j));
                        dRij_detaik          = squeeze(dRi_detai(:,:,k,j));
                        d2epsij_detaikdetail = squeeze(d2epsi_detaidetai(:,k,l,j));
                        d2Rij_detaikdetail   = squeeze(d2Ri_detaidetai(:,:,k,l,j));
                        
                        terms(j) = 2*depsij_detail'*Rij_*depsij_detaik - 2*epsij'*Rij_*dRij_detail*Rij_*depsij_detaik + ...
                            2*epsij'*Rij_*d2epsij_detaikdetail - epsij'*Rij_*d2Rij_detaikdetail*Rij_*epsij + ...
                            2*epsij'*Rij_*dRij_detaik*Rij_*dRij_detail*Rij_*epsij - ...
                            2*epsij'*Rij_*dRij_detaik*Rij_*depsij_detail - ...
                            trace(Rij_*dRij_detail*Rij_*dRij_detaik) + ...
                            trace(Rij_*d2Rij_detaikdetail);
                    end
                    
                    d2li_detai2(k,l) = -1/2 * sum(terms) - Omega_(k,l);
                end
            end
            
            % Eq 47 eta deriv selection term
            detai_detai = eye(neta); % TODO: check this
            
            % Eq 47
            d2li_detaidtheta = zeros(neta,ntheta);
            for k = 1:neta
                for m = 1:ntheta
                    dOmega_dthetam = squeeze(dOmega_dtheta(:,:,m));
                    
                    terms = zeros(ni,1);
                    for j = 1:ni
                        epsij                 = squeeze(epsi(:,j));
                        Rij_                  = squeeze(Ri_(:,:,j));
                        depsij_dthetam        = squeeze(depsi_dtheta(:,m,j));
                        depsij_detaik         = squeeze(depsi_detai(:,k,j));
                        dRij_dthetam          = squeeze(dRi_dtheta(:,:,m,j));
                        dRij_detaik           = squeeze(dRi_detai(:,:,k,j));
                        d2epsij_detaikdthetam = squeeze(d2epsi_detaidtheta(:,k,m,j));
                        d2Rij_detaikdthetam   = squeeze(d2Ri_detaidtheta(:,:,k,m,j));
                        detai_detaik          = squeeze(detai_detai(:,k));
                        
                        terms(j) = 2*depsij_dthetam'*Rij_*depsij_detaik - 2*epsij'*Rij_*dRij_dthetam*Rij_*depsij_detaik + ...
                            2*epsij'*Rij_*d2epsij_detaikdthetam - epsij'*Rij_*d2Rij_detaikdthetam*Rij_*epsij + ...
                            2*epsij'*Rij_*dRij_detaik*Rij_*dRij_dthetam*Rij_*epsij - ...
                            2*epsij'*Rij_*dRij_detaik*Rij_*depsij_dthetam + ...
                            trace(Rij_*dRij_dthetam*Rij_*dRij_detaik + Rij_*d2Rij_detaikdthetam) - ...
                            etai'*Omega_*dOmega_dthetam*Omega_*detai_detaik;
                    end
                    
                    d2li_detaidtheta(k,m) = -1/2 * sum(terms) - etai'*dOmega_dthetam*Omega_*detai_detaik;
                end
            end
            
            % Eq 46 for Eq 33 and Eq 35
            detaistar_dtheta = -inv(d2li_detai2)*d2li_detaidtheta; % neta x ntheta matrix, no time dependence
            
            % Eq 46 corresponding to eta for Eq 35 corresponding to Rij
            detaistar_deta = eye(neta); % neta x neta matrix, no time dependence, selects eta
            
            depsistar_dtheta = zeros(nh,ntheta,ni);
            dRistar_dtheta   = zeros(nh,nh,ntheta,ni);
            dRistar_detai    = zeros(nh,nh,neta,ni);
            for j = 1:ni
                
                dRij_detai   = reshape(dRi_detai(:,:,:,j), nh,nh,neta);
                
                for m = 1:ntheta
                    
                    detaistar_dthetam = squeeze(detaistar_dtheta(:,m));
                    
                    % Eq 33
                    depsij_dthetam    = squeeze(depsi_dtheta(:,m,j));
                    depsij_detai      = squeeze(depsi_detai(:,:,j));
                    depsistar_dtheta(:,m,j) = depsij_dthetam + depsij_detai*detaistar_dthetam;
                    
                    % Eq 35
                    dRij_dthetam = squeeze(dRi_dtheta(:,:,m,j));
                    dRistar_dtheta(:,:,m,j) = dRij_dthetam + tprod(dRij_detai, [1,2,-1], detaistar_dthetam, [-1]);
                end
                
                for k = 1:neta
                    
                    % Eq 35 corresponding to Rij
                    dRij_detaik     = squeeze(dRi_detai(:,:,k,j));
                    detaistar_detak = squeeze(detaistar_deta(:,k));
                    dRistar_detai(:,:,k,j) = dRij_detaik + tprod(dRij_detai, [1,2,-1], detaistar_detak, [-1]);
                end
            end
            
            d_dtheta_depsistar_detai = zeros(nh,ntheta,neta,ni);
            d_dtheta_dRistar_detai   = zeros(nh,nh,ntheta,neta,ni);
            for j = 1:ni
                for m = 1:ntheta
                    for k = 1:neta
                        
                        pdh_dxij           = squeeze(pdh_dxi(:,:,j));
                        pdRij_dxij         = reshape(pdRi_dxi(:,:,:,j), nh,nh,nx);
                        pd2h_detaikdxij    = reshape(pd2h_detaidxi(:,k,:,j), nh,nx);
                        pd2Rij_detaikdxij  = reshape(pd2Ri_detaidxi(:,:,k,:,j), nh,nh,nx);
                        pd2h_dxij2         = reshape(pd2h_dxi2(:,:,:,j), nh,nx,nx);
                        pd2Rij_dxij2       = reshape(pd2Ri_dxi2(:,:,:,:,j), nh,nh,nx,nx);
                        dxij_detaik        = squeeze(dxi_detai(:,k,j));
                        
                        eq37terms = zeros(nh,neta);
                        eq37Rijterms = zeros(nh,nh,neta);
                        for l = 1:neta
                            
                            dxij_detail        = squeeze(dxi_detai(:,l,j));
                            d2xij_detaikdetail = squeeze(d2xi_detaidetai(:,k,l,j));
                            detailstar_dthetam = squeeze(detaistar_dtheta(l,m));
                            
                            pd2h_detaikdetail  = squeeze(pd2h_detaidetai(:,k,l,j));
                            pd2h_dxijdetail    = squeeze(pd2h_dxidetai(:,:,l,j));
                            
                            eq37terms(:,l) = (pd2h_detaikdetail + pd2h_detaikdxij*dxij_detail + ...
                                (pd2h_dxijdetail + tprod(pd2h_dxij2, [1,2,-1], dxij_detail, [-1]))*dxij_detaik + ...
                                pdh_dxij*d2xij_detaikdetail)*detailstar_dthetam;
                            
                            pd2Rij_detaikdetail = squeeze(pd2Ri_detaidetai(:,:,k,l,j));
                            pd2Rij_dxijdetail   = reshape(pd2Ri_dxidetai(:,:,:,l,j), nh,nh,nx);
                            
                            eq37Rijterms(:,:,l) = (pd2Rij_detaikdetail + tprod(pd2Rij_detaikdxij, [1,2,-1], dxij_detail, [-1]) + ...
                                tprod((pd2Rij_dxijdetail + tprod(pd2Rij_dxij2, [1,2,3,-1], dxij_detail, [-1])), [1,2,-1], dxij_detaik, [-1]) + ...
                                tprod(pdRij_dxij, [1,2,-1], d2xij_detaikdetail, [-1]))*detailstar_dthetam;
                        end
                        
                        % Eq 37
                        d_dtheta_depsistar_detai(:,m,k,j) = squeeze(d2epsi_detaidtheta(:,k,m,j)) + sum(eq37terms,2); % first term already has the negative sign
                        
                        % Eq 37 corresponding Rij
                        d_dtheta_dRistar_detai(:,:,m,k,j) = squeeze(d2Ri_detaidtheta(:,:,k,m,j)) + sum(eq37Rijterms,3);
                    end
                end
            end
            
            dastar_dtheta = zeros(nh,neta,ntheta,ni);
            dBstar_dtheta = zeros(nh,nh,ntheta,ni);
            dcstar_dtheta = zeros(nh,nh,neta,ntheta,ni);
            for j = 1:ni
                
                Rijstar_ = squeeze(Ristar_(:,:,j));
                
                for m = 1:ntheta
                    
                    dRijstar_dthetam = squeeze(dRistar_dtheta(:,:,m,j));
                    
                    for k = 1:neta
                        
                        d_dthetam_depsijstar_detaik = squeeze(d_dtheta_depsistar_detai(:,m,k,j));
                        d_dthetam_dRijstar_detaik   = squeeze(d_dtheta_dRistar_detai(:,:,m,k,j));
                        depsijstar_dthetam          = squeeze(depsistar_dtheta(:,m,j));
                        dRijstar_detaik             = squeeze(dRistar_detai(:,:,k,j));
                        epsijstar                   = squeeze(epsistar(:,j));
                        
                        % Eq 30
                        dastar_dtheta(:,k,m,j) = d_dthetam_depsijstar_detaik' - depsijstar_dthetam'*Rijstar_*dRijstar_detaik + ...
                            epsijstar'*Rijstar_*dRijstar_dthetam*Rijstar_*dRijstar_detaik - ...
                            epsijstar'*Rijstar_*d_dthetam_dRijstar_detaik;
                        
                        % Eq 32
                        dcstar_dtheta(:,:,k,m,j) = -Rijstar_*dRijstar_dthetam*Rijstar_*dRijstar_detaik + ...
                            Rijstar_*d_dthetam_dRijstar_detaik;
                    end
                    
                    % Eq 31
                    dBstar_dtheta(:,:,m,j) = -2*Rijstar_*dRijstar_dthetam*Rijstar_;
                end
            end
            
            % Eq 28
            dli_dtheta = zeros(ntheta,1);
            for m = 1:ntheta
                
                eq28terms = zeros(ni,1); % scalar
                for j = 1:ni
                    epsijstar          = squeeze(epsistar(:,j));
                    Rijstar_           = squeeze(Ristar_(:,:,j));
                    depsijstar_dthetam = squeeze(depsistar_dtheta(:,m,j));
                    dRijstar_dthetam   = squeeze(dRistar_dtheta(:,:,m,j));
                    
                    eq28terms(j) = 2*epsijstar'*Rijstar_*depsijstar_dthetam - ...
                        epsijstar'*Rijstar_*dRijstar_dthetam*Rijstar_*epsijstar + ...
                        trace(Rijstar_*dRijstar_dthetam);
                end
                
                dOmega_dthetam = squeeze(dOmega_dtheta(:,:,m));
                
                dli_dtheta(m) = -1/2 * sum(eq28terms) + ...
                    1/2 * etai'*Omega_*dOmega_dthetam*Omega_*etai - ...
                    1/2 * trace(Omega_*dOmega_dthetam);
            end
            
            % Eq 29
            dHi_dtheta = zeros(neta,neta,ntheta);
            for k = 1:neta
                for l = 1:neta
                    for m = 1:ntheta
                        dHikl = zeros(ni,1); % term inside the sum
                        for j = 1:ni
                            dalstar_dthetamj = squeeze(dastar_dtheta(:,l,m,j));
                            dakstar_dthetamj = squeeze(dastar_dtheta(:,k,m,j));
                            Bstarj = squeeze(Bstar(:,:,j));
                            dBstar_dthetamj = squeeze(dBstar_dtheta(:,:,m,j));
                            akstarj = squeeze(astar(:,k,j));
                            alstarj = squeeze(astar(:,l,j));
                            dclstar_dthetamj = squeeze(dcstar_dtheta(:,:,l,m,j));
                            dckstar_dthetamj = squeeze(dcstar_dtheta(:,:,k,m,j));
                            ckstarj = squeeze(cstar(:,:,k,j));
                            clstarj = squeeze(cstar(:,:,l,j));
                            
                            dHikl(j) = dalstar_dthetamj*Bstarj*akstarj' + ...
                                alstarj*dBstar_dthetamj*akstarj' + ...
                                alstarj*Bstarj*dakstar_dthetamj' + ...
                                trace(-dclstar_dthetamj*ckstarj - clstarj*dckstar_dthetamj);
                        end
                        
                        dOmegakl_dthetam_ = squeeze(dOmega_dtheta_(k,l,m));
                        
                        dHi_dtheta(k,l,m) = -1/2 * sum(dHikl) - dOmegakl_dthetam_;
                    end
                end
            end
            
            % Eq 23, for each patient
            Hi_ = inv(Hi);
            dlogLFi_dthetam = zeros(ntheta,1);
            for m = 1:ntheta
                dHi_dthetam = squeeze(dHi_dtheta(:,:,m));
                dlogLFi_dthetam(m) = dli_dtheta(m) - 1/2 * trace(Hi_*dHi_dthetam);
            end
            
            % Assign output obj fun gradient
            dGdk = dlogLFi_dthetam
            
        end
        
    end

    function [val, discrete] = G(int)
        order = 1;
        [val, discrete] = calcParts(int, order); %%%DEBUG%%% look at curvature for now
        %         val = 0;
        %         discrete = 0;
    end

    function val = dGdk(int)
        % This is called when calculating the obj fun val only and isn't needed.
        % When it is needed, as when calculating the obj fun grad, int will be
        % augmented with the 2nd order derivatives.
        order = 2;
        if ~isfield(int, 'd2ydx2')
            val = zeros(int.nk, 1);
        else
            val = calcParts(int, order);
        end
    end

    function val = d2Gdk2(int)
        % This is called when calculating the obj gradient and isn't needed
        %warning('Analytic derivative for NLME d2Gdk2 not implemented yet. Derive and implement it or approximate with finite differences.')
        order = 3;
        if ~isfield(int, 'd3ydx3') % will never be called
            val = zeros(int.nk, int.nk);
        else
            val = calcParts(int, order);
        end
    end


end

%% Helper functions
function [y_Ind, Fi_y_Ind] = lookupOutput(name, int)
y_Ind    = find(ismember(int.y_names, name));
Fi_y_Ind = find(ismember(int.y_names, ['Fi__' name]));
end

function [RiInds, hInds, errInds] = getOutputInds(yNames)
ny = length(yNames);
RiInds  = zeros(ny,1);
hInds   = zeros(ny,1);
errInds = zeros(ny,1);
for i = 1:ny
    yName = yNames{i};
    if ~isempty(regexp(yName, '^Ri__', 'once'))
        RiInds(i) = i;
    elseif ~isempty(regexp(yName, '^h__', 'once'))
        hInds(i) = i;
    elseif ~isempty(regexp(yName, '^err__', 'once'))
        errInds(i) = i;
    else
        % not recognized, ignore
    end
end
RiInds(RiInds==0)   = [];
hInds(hInds==0)     = [];
errInds(errInds==0) = [];
end

function [thetaInds, etaInds, omegaInds, sigmaInds] = getNLMEInds(kNames)
% Get indices of NLME components from rate params list
nk = length(kNames);
thetaInds = zeros(nk,1);
etaInds   = zeros(nk,1);
omegaInds = zeros(nk,1);
sigmaInds = zeros(nk,1);
for i = 1:nk
    kName = kNames{i};
    if ~isempty(regexp(kName, '^eta__', 'once'))
        etaInds(i) = i;
    elseif ~isempty(regexp(kName, '^omega__', 'once'))
        omegaInds(i) = i;
    elseif ~isempty(regexp(kName, '^sigma__', 'once'))
        sigmaInds(i) = i;
    else % not a NLME param, so call a theta
        thetaInds(i) = i;
    end
end
thetaInds(thetaInds==0) = [];
etaInds(etaInds==0)     = [];
omegaInds(omegaInds==0) = [];
sigmaInds(sigmaInds==0) = [];
end

function Omega = getOmega(k, kNames, thetaInds, omegaInds)
% Assemble full Omega matrix by looking up params from k and kNames. Omega is a
% symmetric n x n matrix where n = number of "eta" parameters being fit and
% order the same as k.
% Parameters w/o inter-individual variability eta specified have omegas set = 0.

omegas = k(omegaInds);
omegaNames = kNames(omegaInds);
thetaNames = kNames(thetaInds);

no = length(omegaInds);
n = length(thetaInds);
Omega = zeros(n);
for i = 1:no
    tokens = regexp(omegaNames{i}, '^omega__(.+)__(.+)$', 'tokens', 'once')';
    assert(iscell(tokens) && length(tokens) == 2, 'observationNLME:getOmega:invalidOmegaTokens', 'Parser didn''t find correct params omega name %s is referring to', omegaNames{i})
    ix = find(ismember(thetaNames, tokens{1}));
    jx = find(ismember(thetaNames, tokens{2}));
    
    Omega(ix,jx) = omegas(i);
    Omega(jx,ix) = omegas(i);
end
end

% function Sigma = getSigma(k, kNames, epsInds, sigmaInds)
% % Assemble full Sigma matrix by looking up params from k and kNames. Sigma is a
% % symmetric n x n matrix where n = number of "eps" parameters being fit.
%
% sigmas = k(sigmaInds);
% sigmaNames = kNames(sigmaInds);
% epsNames   = kNames(epsInds);
%
% % Trim leading "eps__" from epsNames
% n = length(epsInds);
% for i = 1:n
%     epsNames{i} = strrep(epsNames{i}, 'eps__', '');
% end
%
% no = length(sigmaInds);
% Sigma = zeros(n);
% for i = 1:no
%     tokens = regexp(sigmaNames{i}, '^sigma__(.+)__(.+)$', 'tokens', 'once')';
%     assert(iscell(tokens) && length(tokens) == 2, 'observationNLME:getSigma:invalidSigmaTokens', 'Parser didn''t find correct params sigma name %s is referring to', sigmaNames{i})
%     ix = find(ismember(epsNames, tokens{1}));
%     jx = find(ismember(epsNames, tokens{2}));
%
%     Sigma(ix,jx) = sigmas(i);
%     Sigma(jx,ix) = sigmas(i);
% end
% end

% function Gi = getGi(dydT, kNames, thetaInds, etaInds, Fi_y_Ind, ny)
% % Gi = dF/deta, a nTimes (nObservations) x netas matrix
% % dydT represents the ny x nT x nt matrix of partial derivatives of
% %   outputs wrt all parameters being optimized over as (ny x nT) * nt. Going
% %   down dydT, y changes fastest, i.e., the first ny entries are dy1/dT1, dy2/dT1, dy3/dT1, etc.
% %   the next ny entries are dy1/dT2, etc.
% % kNames is sufficient because k's are first in T and other types of params aren't used
%
% thetaNames = kNames(thetaInds);
% etaNames   = kNames(etaInds);
%
% nTheta = length(thetaInds);
% nEta = length(etaInds);
%
% % Get ks corresponding to etas
% dydT_eta_Ind = zeros(nTheta,1);
% for i = 1:nEta
%     token = regexp(etaNames{i}, '^eta__(.+)$', 'tokens', 'once');
%     thetaNameThisEta = token{1};
%     thetaMaskThisEta = ismember(thetaNames, thetaNameThisEta);
%     dydT_eta_Ind(thetaMaskThisEta) = find(ismember(kNames, etaNames{i}));
% end
%
% % Assemble Gi matrix
% % Thetas with no corresponding etas will have 0 values
% nt = size(dydT,2);
% Gi = zeros(nt,nTheta);
% for i = 1:nTheta
%     offset = ny*(dydT_eta_Ind(i)-1); % to find the right parameter block
%     ind = Fi_y_Ind + offset; % to find the right output
%     if ind ~= 0
%         Gi(:,i) = dydT(ind,:)';
%     end
% end
% end

function Hi = getHi(dydT, epsInds, y_Ind, ny)
% Hi = dY/deps
nEps = length(epsInds);
nt = size(dydT,2);
Hi = zeros(nt,nEps);
for i = 1:nEps
    offset = ny*(epsInds(i)-1); % to find the right parameter block
    ind = y_Ind + offset; % to find the right output
    Hi(:,i) = dydT(ind,:)';
end
end

function output = fixOutputName(output)
% Convert output argument to right form
if iscell(output)
    % pass
elseif ischar(output)
    output = {output};
else
    error('observationNLME:invalidOutputArg', 'Output must be a string or cell array of strings')
end
nOutputs = length(output);
% Prefix outputs with 'Ri__' if not already
for i = 1:nOutputs
    if ~ischar(output{i})
        error('observationNLME:invalidOutputType', 'Output must be a cell array of strings')
    end
    if ~strncmp(output{i}, 'Ri__', 4)
        output{i} = ['Ri__' output{i}];
    end
end
end