function [objfunval, objfungrad] = calcFocei(int, outputs, measurements, order, layer)
% Calculate components of NLME obj fun
%
% Inputs:
%   int [ struct ]
%       Integrated solution struct
%   outputs [ nOutputs x 1 cell vector of strings ]
%       Names of outputs
%   measurements [ nOutputs x ni double matrix ]
%       Individuals' measurements for each output at each time
%   order [ 1 | 2 ]
%       Whether to require sensitivity (1) or curvature (2) to
%       calculate obj fun value or gradient. Redundant
%       check with presence of int.d2ydx2
%   layer [ 'outer' | 'inner' ]
%       Whether to return outer optimization obj fun val and
%       gradient (logli and dlogli_dtheta) or inner optimization obj fun
%       val and gradient (li and dli_deta)
%
% Outputs:
%   objfunval [ scalar double ]
%       
%   objfungrad [ ntheta x 1 double vector or neta x 1 double vector ]
%       
%
% Notes:
%   Only use squeeze when final result is 1 or 2 dimensional and kept
%   dimensions are adjacent
%       Using reshape can remove singleton dimensions unnecessarily when
%   Careful with squeezing/reshaping Ri with single output when tensor is needed - will remove extra singleton dims
%       Only call squeeze when colons are adjacent and in order, else call reshape
%   Partial derivatives from function handles; total derivative from filled in values
%       partial derivatives start with pd
%   Fixed effects parameters = Theta, random effects parameters = Eta
%       Omega and Sigma are included in Theta
%       Omega is pulled out separately for convenience for a couple terms
%   All terms with star are evaluated at optimal eta star (same as eta in FO method)
%   int.dxdT and int.dx2dT2 are based on T, the vector of params this condition
%       is fitting (and not all of them). For an individual, use UseParams and
%       UseSeeds (TODO: implement UseInputControls and UseDoseControls also) to
%       do indexing correctly

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get common components
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nx = int.nx;
ny = int.ny;
nT = int.nT;
% if isfield(int, 'dxdT')
%     nT = int.nT;
% end

% Indexing
[hInds, RiInds, RiPos] = getOutputInds(int.y_names, outputs);
nh = length(hInds);
nRi = length(RiInds);

[thetaInds, etaInds] = getParameterInds(int.k_names);
ntheta = length(thetaInds);
neta = length(etaInds);

% Indices that take into account the subset of fit params when using T
[thetaIndsT, etaIndsT] = getParameterInds(int.k_names(int.UseParams));
nthetaT = length(thetaIndsT);
netaT = length(etaIndsT);

etai = int.k(etaInds);

% Do everything by timepoints
ni = size(measurements,2);

if order >= 1
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Assemble obj fun val components
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    epsi        = zeros(nh,ni);
    
    Ri          = zeros(nh,nh,ni);
    Ri_         = zeros(nh,nh,ni);
    
    pdh_detai   = zeros(nh,neta,ni);
    pdh_dxi     = zeros(nh,nx,ni);
    pdRi_detai  = zeros(nh,nh,neta,ni);
    pdRi_dxi    = zeros(nh,nh,nx,ni);
    dxi_detai   = zeros(nx,neta,ni);
    
    Omega       = int.fOmega.Omega(int.k);
    Omega_      = int.fOmega.OmegaI(int.k);
    
    for j = 1:ni
        
        % Eq 5
        %   Measurements is a ni x nh matrix (transposed from everything else)
        %   yhat used here is E[y], excludes intra-individual error term
        epsij = zeros(nh,1);
        for ih = 1:nh
            epsij(ih) = measurements(ih,j) - int.y(hInds(ih),j);
        end
        epsi(:,j) = epsij;
        
        % Partial derivatives from function handles, time dependence from
        %   specifying x and u at time indices
        pdydxj = int.dydx([], int.x(:,j), int.u(:,j));
        pdydkj = int.dydk([], int.x(:,j), int.u(:,j));
        
        % Get Ri, intra-individual variability model covariance matrix
        % Get inverse of Ri once for convenience
        Rij = zeros(nh);
        for iRi = 1:nRi
            ih = RiPos(iRi,1);
            jh = RiPos(iRi,2);
            Rij(ih,jh) = int.y(RiInds(iRi),j);
            Rij(jh,ih) = int.y(RiInds(iRi),j);
        end
        Rij_ = inv(Rij);
        Ri(:,:,j)  = Rij;
        Ri_(:,:,j) = Rij_;
        
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
        for iRi = 1:nRi
            ih = RiPos(iRi,1);
            jh = RiPos(iRi,2);
            for k = 1:neta
                pdRij_detaik(ih,jh,k) = pdydkj(RiInds(iRi),etaInds(k));
                pdRij_detaik(jh,ih,k) = pdydkj(RiInds(iRi),etaInds(k));
            end
        end
        pdRi_detai(:,:,:,j) = pdRij_detaik;
        
        % In Eq 20
        pdRij_dxij = zeros(nh,nh,nx);
        for iRi = 1:nRi
            ih = RiPos(iRi,1);
            jh = RiPos(iRi,2);
            for ix = 1:nx
                pdRij_dxij(ih,jh,ix) = pdydxj(RiInds(iRi),ix);
                pdRij_dxij(jh,ih,ix) = pdydxj(RiInds(iRi),ix);
            end
        end
        pdRi_dxi(:,:,:,j) = pdRij_dxij;
        
        % In Eq 19,20
        dxij_detaik = zeros(nx,neta);
        for ix = 1:nx
            for k = 1:netaT
                ind = nx*(etaIndsT(k)-1) + ix; %%FIT_T_ONLY%% what are the etaInds in _dt
                dxij_detaik(ix,k) = int.dxdT(ind,j);
            end
        end
        dxi_detai(:,:,j) = dxij_detaik;
        
    end
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate obj fun val
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    a           = zeros(nh,neta,ni);
    B           = zeros(nh,nh,ni);
    c           = zeros(nh,nh,neta,ni);
    
    depsi_detai = zeros(nh,neta,ni);
    dRi_detai   = zeros(nh,nh,neta,ni);
    
    for j = 1:ni
        for k = 1:neta
            
            % Eq 19
            pdh_detaik  = squeeze(pdh_detai(:,k,j));
            pdh_dxij    = squeeze(pdh_dxi(:,:,j));
            dxij_detaik = squeeze(dxi_detai(:,k,j));
            depsij_detaik = -(pdh_detaik + pdh_dxij*dxij_detaik);
            depsi_detai(:,k,j) = depsij_detaik;
            
            % Eq 20
            pdRij_detaik = squeeze(pdRi_detai(:,:,k,j));
            pdRij_dxij   = reshape(pdRi_dxi(:,:,:,j), nh,nh,nx);
            dRij_detaik = pdRij_detaik + tprod(pdRij_dxij, [1,2,-1], dxij_detaik, [-1]); % contract on xi
            dRi_detai(:,:,k,j) = dRij_detaik;
            
            % Eq 15
            epsij = squeeze(epsi(:,j));
            Rij_  = squeeze(Ri_(:,:,j));
            a(:,k,j) = depsij_detaik' - epsij'*Rij_*dRij_detaik;
            
            % Eq 16
            B(:,:,j) = 2*Rij_;
            
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
                al = reshape(a(:,l,j), 1,nh);
                ak = reshape(a(:,k,j), 1,nh);
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
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Select output for outer or inner optimization
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(layer, 'outer')
        % Eq 18, for each patient
        logLFi = li - 1/2 * log(det(-Hi/(2*pi))); % Almquist
        %             logLFi = li - 1/2 * log(det(-Hi_tilde/(2*pi))); % NONMEM
        objfunval = -logLFi; % negate to maximize log-likelihood
    elseif strcmp(layer, 'inner')
        objfunval = -li; % negate to maximize likelihood
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obj fun grad calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if order >= 2
    
    if strcmp(layer, 'inner')
        % Inner optimization obj fun grad only require 1st order
        %   sensitivities
        
        % Eq 12
        dli_detai = zeros(neta,1);
        for k = 1:neta
            
            eq12terms = zeros(ni,1);
            for j = 1:ni
                epsij = squeeze(epsi(:,j));
                Rij_  = squeeze(Ri_(:,:,j));
                depsij_detaik = squeeze(depsi_detai(:,k,j));
                dRij_detaik = squeeze(dRi_detai(:,:,k,j));
                
                eq12terms(j) = 2*epsij'*Rij_*depsij_detaik - ...
                    epsij'*Rij_*dRij_detaik*Rij_*epsij + ...
                    trace(Rij_*dRij_detaik);
            end
            
            detai_detai = zeros(neta,1);
            detai_detai(k) = 1;
            dli_detai(k) = -1/2*sum(eq12terms) - etai'*Omega_*detai_detai;
        end
        objfungrad = -dli_detai; % negate to maximize likelihood
        
    elseif strcmp(layer, 'outer')
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Assemble obj fun grad components
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Assume/TODO - these may be just the correct values when int is reevaluated at etastar
        Ristar  = Ri;
        Ristar_ = Ri_;
        
        dOmega_dtheta  = int.fOmega.dOmega_dtheta(int.k);
        dOmega_dtheta_ = int.fOmega.dOmegaI_dtheta(int.k);
        dOmega_dtheta(:,:,etaInds)  = []; % remember to remove thetas corresponding to etas
        dOmega_dtheta_(:,:,etaInds) = [];
        
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
            for iRi = 1:nRi
                ih = RiPos(iRi,1);
                jh = RiPos(iRi,2);
                for m = 1:ntheta
                    pdRij_dtheta(ih,jh,m) = pdydkj(RiInds(iRi),thetaInds(m));
                    pdRij_dtheta(jh,ih,m) = pdydkj(RiInds(iRi),thetaInds(m));
                end
            end
            pdRi_dtheta(:,:,:,j) = pdRij_dtheta;
            
            % In Eq 34
            dxij_dtheta = zeros(nx,ntheta);
            for ix = 1:nx
                for m = 1:nthetaT
                    ind = nx*(thetaIndsT(m)-1) + ix;
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
            for ih = 1:nh
                for k = 1:neta
                    for m = 1:ntheta
                        etaInd = etaInds(k);
                        hthetaInd = ny*(thetaInds(m)-1) + hInds(ih);
                        pd2h_detaidtheta(ih,k,m) = pd2ydk2j(hthetaInd,etaInd);
                    end
                end
            end
            pd2h_detaidtheta(:,:,:,j) = pd2h_detaidtheta;
            
            pd2Rij_detaidtheta = zeros(nh,nh,neta,ntheta);
            for iRi = 1:nRi
                ih = RiPos(iRi,1);
                jh = RiPos(iRi,2);
                for k = 1:neta
                    for m = 1:ntheta
                        etaInd = etaInds(k);
                        RithetaInd = ny*(thetaInds(m)-1) + RiInds(iRi);
                        
                        pd2Rij_detaidtheta(ih,jh,k,m) = pd2ydk2j(RithetaInd,etaInd);
                        pd2Rij_detaidtheta(jh,ih,k,m) = pd2ydk2j(RithetaInd,etaInd);
                    end
                end
            end
            pd2Ri_detaidtheta(:,:,:,:,j) = pd2Rij_detaidtheta;
            
            % In Eq 37
            pd2h_detaidxij = zeros(nh,neta,nx);
            for ih = 1:nh
                for k = 1:neta
                    for ix = 1:nx
                        etaInd = etaInds(k);
                        hxInd = ny*(ix-1) + hInds(ih);
                        pd2h_detaidxij(ih,k,ix) = pd2ydkdxj(hxInd,etaInd);
                    end
                end
            end
            pd2h_detaidxi(:,:,:,j) = pd2h_detaidxij;
            
            pd2Rij_detaidxij = zeros(nh,nh,neta,nx);
            for iRi = 1:nRi
                ih = RiPos(iRi,1);
                jh = RiPos(iRi,2);
                for k = 1:neta
                    for ix = 1:nx
                        etaInd = etaInds(k);
                        RixInd = ny*(ix-1) + RiInds(iRi);
                        pd2Rij_detaidxij(ih,jh,k,ix) = pd2ydkdxj(RixInd,etaInd);
                        pd2Rij_detaidxij(jh,ih,k,ix) = pd2ydkdxj(RixInd,etaInd);
                    end
                end
            end
            pd2Ri_detaidxi(:,:,:,:,j) = pd2Rij_detaidxij;
            
            % In Eq 37
            pd2h_dxijdtheta = zeros(nh,nx,ntheta);
            for ih = 1:nh
                for ix = 1:nx
                    for m = 1:ntheta
                        xInd = ix;
                        hthetaInd = ny*(thetaInds(m)-1) + hInds(ih);
                        pd2h_dxijdtheta(ih,ix,m) = pd2ydxdkj(hthetaInd,xInd);
                    end
                end
            end
            pd2h_dxidtheta(:,:,:,j) = pd2h_dxijdtheta;
            
            pd2Rij_dxijdtheta = zeros(nh,nh,nx,ntheta);
            for iRi = 1:nRi
                ih = RiPos(iRi,1);
                jh = RiPos(iRi,2);
                for ix = 1:nx
                    for m = 1:ntheta
                        xInd = ix;
                        RithetaInd = ny*(thetaInds(m)-1) + RiInds(iRi);
                        pd2Rij_dxijdtheta(ih,jh,ix,m) = pd2ydxdkj(RithetaInd,xInd);
                        pd2Rij_dxijdtheta(jh,ih,ix,m) = pd2ydxdkj(RithetaInd,xInd);
                    end
                end
            end
            pd2Ri_dxidtheta(:,:,:,:,j) = pd2Rij_dxijdtheta;
            
            % In Eq 37
            pd2h_dxij2 = zeros(nh,nx,nx);
            for ih = 1:nh
                for ix = 1:nx
                    for jx = 1:nx
                        xInd = ix;
                        hxInd = ny*(jx-1) + hInds(ih);
                        pd2h_dxij2(ih,ix,jx) = pd2ydx2j(hxInd,xInd);
                    end
                end
            end
            pd2h_dxi2(:,:,:,j) = pd2h_dxij2;
            
            pd2Rij_dxij2 = zeros(nh,nh,nx,nx);
            for iRi = 1:nRi
                ih = RiPos(iRi,1);
                jh = RiPos(iRi,2);
                for ix = 1:nx
                    for jx = 1:nx
                        xInd = ix;
                        RixInd = ny*(jx-1) + RiInds(iRi);
                        pd2Rij_dxij2(ih,jh,ix,jx) = pd2ydx2j(RixInd,xInd);
                        pd2Rij_dxij2(jh,ih,ix,jx) = pd2ydx2j(RixInd,xInd);
                    end
                end
            end
            pd2Ri_dxi2(:,:,:,:,j) = pd2Rij_dxij2;
            
            % In Eq 37
            pd2h_detaidetaij = zeros(nh,neta,neta);
            for ih = 1:nh
                for k = 1:neta
                    for l = 1:neta
                        etaInd = etaInds(k);
                        hetaInd = ny*(etaInds(l)-1) + hInds(ih);
                        pd2h_detaidetaij(ih,k,l) = pd2ydk2j(hetaInd,etaInd);
                    end
                end
            end
            pd2h_detaidetai(:,:,:,j) = pd2h_detaidetaij;
            
            pd2Rij_detaidetai = zeros(nh,nh,neta,neta);
            for iRi = 1:nRi
                ih = RiPos(iRi,1);
                jh = RiPos(iRi,2);
                for k = 1:neta
                    for l = 1:neta
                        etaInd = etaInds(k);
                        RietaInd = ny*(etaInds(l)-1) + RiInds(iRi);
                        pd2Rij_detaidetai(ih,jh,k,l) = pd2ydk2j(RietaInd,etaInd);
                        pd2Rij_detaidetai(jh,ih,k,l) = pd2ydk2j(RietaInd,etaInd);
                    end
                end
            end
            pd2Ri_detaidetai(:,:,:,:,j) = pd2Rij_detaidetai;
            
            % In Eq 37
            pd2h_dxijdetai = zeros(nh,nx,neta);
            for ih = 1:nh
                for ix = 1:nx
                    for l = 1:neta
                        xInd = ix;
                        hetaInd = ny*(etaInds(l)-1) + hInds(ih);
                        pd2h_dxijdetai(ih,ix,l) = pd2ydxdkj(hetaInd,xInd);
                    end
                end
            end
            pd2h_dxidetai(:,:,:,j) = pd2h_dxijdetai;
            
            pd2Rij_dxijdetai = zeros(nh,nh,nx,neta);
            for iRi = 1:nRi
                ih = RiPos(iRi,1);
                jh = RiPos(iRi,2);
                for ix = 1:nx
                    for l = 1:neta
                        xInd = ix;
                        RietaInd = ny*(etaInds(l)-1) + RiInds(iRi);
                        pd2Rij_dxijdetai(ih,jh,ix,l) = pd2ydxdkj(RietaInd,xInd);
                        pd2Rij_dxijdetai(jh,ih,ix,l) = pd2ydxdkj(RietaInd,xInd);
                    end
                end
            end
            pd2Ri_dxidetai(:,:,:,:,j) = pd2Rij_dxijdetai;
            
            % In Eq 37
            % x changes fastest, then thetam; etaik changes slowest
            d2xij_detaikdthetam = zeros(nx,neta,ntheta);
            for ix = 1:nx
                for k = 1:netaT
                    for m = 1:nthetaT
                        ind = nx*nT*(etaIndsT(k)-1) + nx*(thetaIndsT(m)-1) + ix;
                        d2xij_detaikdthetam(ix,k,m) = int.d2xdT2(ind,j); % 2 x 13 x 13 x time; 1st 3 dims compressed into row indices
                    end
                end
            end
            d2xi_detaidtheta(:,:,:,j) = d2xij_detaikdthetam;
            
            % In Eq 37
            % x changes fastest, then etail; etaik changes slowest
            d2xij_detaikdetail = zeros(nx,neta,neta);
            for ix = 1:nx
                for k = 1:netaT
                    for l = 1:netaT
                        ind = nx*nT*(etaIndsT(k)-1) + nx*(etaIndsT(l)-1) + ix;
                        d2xij_detaikdetail(ix,k,l) = int.d2xdT2(ind,j); % 2 x 13 x 13 x time; 1st 3 dims compressed into row indices
                    end
                end
            end
            d2xi_detaidetai(:,:,:,j) = d2xij_detaikdetail;
        end
        
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Calculate obj fun grad
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Assume/TODO - these may be just the correct values when int is reevaluated at etastar
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
                    pd2Rij_dxijdthetam   = reshape(pd2Ri_dxidtheta(:,:,:,m,j), nh,nh,nx);
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
        detaistar_deta = eye(neta); % neta x neta matrix, no time dependence, selects eta - TODO/check this
        
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
                %                     depsijstar_dthetam = squeeze(depsistar_dtheta(:,m,j));
                depsijstar_dthetam = squeeze(depsi_dtheta(:,m,j)); % FIX - the term w/o the star is the right one
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
                        dalstar_dthetamj = reshape(dastar_dtheta(:,l,m,j), 1,nh);
                        dakstar_dthetamj = reshape(dastar_dtheta(:,k,m,j), 1,nh);
                        Bstarj           = squeeze(Bstar(:,:,j));
                        dBstar_dthetamj  = squeeze(dBstar_dtheta(:,:,m,j));
                        akstarj          = reshape(astar(:,k,j), 1,nh);
                        alstarj          = reshape(astar(:,l,j), 1,nh);
                        dclstar_dthetamj = squeeze(dcstar_dtheta(:,:,l,m,j));
                        dckstar_dthetamj = squeeze(dcstar_dtheta(:,:,k,m,j));
                        ckstarj          = squeeze(cstar(:,:,k,j));
                        clstarj          = squeeze(cstar(:,:,l,j));
                        
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
        dlogLFi_dtheta = zeros(ntheta,1);
        for m = 1:ntheta
            dHi_dthetam = squeeze(dHi_dtheta(:,:,m));
            dlogLFi_dtheta(m) = dli_dtheta(m) - 1/2 * trace(Hi_*dHi_dthetam);
        end
        
        % Assign output obj fun gradient
        objfungrad = -dlogLFi_dtheta; % negate to maximize likelihood
        
        %% Debug: Obj fun value w/ full Hessian (Eq 13) for comparison
        % Eq 18
        %             logLFiLapacian = li - 1/2 * log(det(-d2li_detai2/(2*pi)));
        %             GLaplacian = -logLFiLapacian;
    end
end
end
