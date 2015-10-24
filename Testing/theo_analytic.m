% Analytic solution and explicit NLME expressions for testing/comparison for
% Theophylline model
% From Min-Gul Kim, et al., 2015
% Note that positions of k = kel and V are switched
function [objValTot, objGradTot] = theo_analytic
theoph = load_theo();

% Same starting params as test case
Theta = [1.5 0.15 25];
Eta = [0 0 0];
Omega = [0.2 0.1 0.1; 0.1 0.2 0.1; 0.1 0.1 0.2];
Sigma = [0.1 0; 0 0.1];

kth = Theta';
kom = Omega(tril(true(size(Omega))));
ksg = Sigma(tril(true(size(Sigma))));
k = [kth; kom; ksg]; % [theta:1:3, omega:4:9, sigma:10:12]

nkth = length(kth);
nkom = length(kom);
nksg = length(ksg);
nk = length(k);

nPatients = height(theoph);

objVal  = zeros(1,nPatients);
objGrad = zeros(nk,nPatients);
for i = 1:nPatients
    data = theoph(i,:);
    objVal(i) = G(k);
%     objGrad(:,i) = dGdk(k);
    objGrad(:,i) = dGdk_finite(k);
end
objValTot  = sum(objVal);
objGradTot = sum(objGrad);

%% Helper functions

    function objVal = obj(k)
        % Calculate objective value
        
        % Extract params
        Theta_i = k(1:nkth);
        Omega_i = ltm2full(ltv2ltm(k(nkth+1:nkth+nkom)));
        Sigma_i = ltm2full(ltv2ltm(k(nkth+nkom+1:nk)));
        
        [Fi, Gi, Hi] = pred(Theta_i, Eta, data);
        Yi = data.Conc{1};
        
        RES = Yi - Fi;
        COV = Gi*Omega_i*Gi' + diag(diag(Hi*Sigma_i*Hi'));
        WRES = matSqrtInv(COV)*RES;
        
        objVal = log(det(COV)) + WRES'*WRES;
        
    end


%% Obj fun and gradient calculations of single patient

    function val = G(k)
        val = obj(k);
    end

    function val = dGdk(k)
        % Order of k is [Theta; Omega; Sigma]
        val = zeros(nk,1);
    end

    function val = dGdk_finite(k)
        % Central finite difference approximation to dGdk
        delta = 1e-6;
        val = zeros(nk,1);
        
        for j = 1:nk
            kj = [k,k];
            kj(j,1) = k(j) - delta/2;
            kj(j,2) = k(j) + delta/2;
            val(j) = (obj(kj(:,2)) - obj(kj(:,1)))/delta;
        end
        
    end

end

%% Helper functions
function [Fi, Gi, Hi] = pred(Theta, Eta, Data)
% Per-patient prediction
dose = 320;
time = Data.Time{1};

ka  = Theta(1)*exp(Eta(1));
kel = Theta(2)*exp(Eta(2));
V   = Theta(3)*exp(Eta(3));

Fi = dose/V * ka/(ka - kel) * (exp(-kel*time) - exp(-ka*time));

G1 = ((dose/V/(ka - kel) - dose/V .* ka/(ka - kel).^2) .* (exp(-kel*time) - exp(-ka*time)) + ...
    dose/V .* ka/(ka - kel) .* (exp(-ka*time).*time)) .* ka;

G2 = (dose/V .* ka/(ka - kel)^2 .* (exp(-kel*time) - exp(-ka*time)) - ...
    dose/V .* ka/(ka - kel) .* (exp(-kel*time).*time)) .* kel; % G3 in paper

G3 = -(dose/V^2 .* ka/(ka - kel) .* (exp(-kel*time) - exp(-ka*time))) .* V; % G2 in paper

H1 = Fi;

H2 = ones(size(Fi));

Gi = [G1, G2, G3];
Hi = [H1, H2];
end

function ltm = ltv2ltm(ltv)
% Convert vector of lower triangular matrix elements to lower triangular matrix
n = round((sqrt(8*length(ltv)+1)-1)/2); % size of output matrix 
ltm = zeros(n);
[ii, jj] = ndgrid(1:n);
ltm(ii>=jj) = ltv; % fill in values in column-major order
end

function full = ltm2full(ltm)
% Convert lower triangular matrix to full symmetric matrix
full = ltm + ltm' - eye(size(ltm)).*ltm;
end