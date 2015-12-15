% Analytic/symbolic version of 2 observation simple model

function T08b2_FO_Linear_2_Analytic
clear; close all; clc
rng('default');

%% Model
% Time
t = sym('t');

% Seeds
x1_0 = sym('x1_0');
x2_0 = sym('x2_0');
x0 = [x1_0 x2_0].';
x0_vals = [10 0.1]';

% Parameters
k1_0 = sym('k1_0');
k2_0 = sym('k2_0');
k_0 = [k1_0 k2_0].';
k_0_vals = [0.9 1.0]';

% Inter-individual parameters
eta_k1 = sym('eta_k1');
eta_k2 = sym('eta_k2');
eta = [eta_k1 eta_k2].';
eta_vals = [0 0]';
neta = length(eta);

% Expressions for parameters with inter-individual variability model
k1 = k1_0*exp(eta_k1);
k2 = k2_0*exp(eta_k2);
k = [k1 k2].';

% States over time (explicit analytic solution)
h1 = x1_0*exp(-k1*t);
h2 = x1_0 - x1_0*exp(-k1*t) + x2_0;
h = [h1 h2]; % row vector
nh = length(h);

% Inter-individual covariance - specified as the omegas^2 variances
omega__k1__k1 = sym('omega__k1__k1');
omega__k2__k1 = sym('omega__k2__k1');
omega__k2__k2 = sym('omega__k2__k2');
Omega = [omega__k1__k1, omega__k2__k1; omega__k2__k1, omega__k2__k2];
omega = [omega__k1__k1 omega__k2__k1 omega__k2__k2].';
omega_vals = [0.1 0.01 0.1]';

% Intra-individual covariance - specified as individual params sigmas
sigma__y1 = sym('sigma__y1');
sigma__y2__y1 = sym('sigma__y2__y1');
sigma__y2 = sym('sigma__y2');
Ri = [sigma__y1.^2, sigma__y2__y1.^2; sigma__y2__y1.^2, sigma__y2.^2];
Ri_ = inv(Ri); % for convenience
sigma = [sigma__y1 sigma__y2__y1 sigma__y2].';
sigma_vals = [0.05 0.01 0.04]';

% Combined parameter vector
theta = [k1_0 k2_0 omega__k1__k1 omega__k2__k1 omega__k2__k2 sigma__y1 sigma__y2__y1 sigma__y2].';
ntheta = length(theta);

%% Simulate some experimental observations
ni = 5; % number of observations/times
t0 = 0;
tf = 10;
times = sym('t', [ni,1]);
tvals = linspace(t0, tf, ni)';

d = sym('d', [ni,2]);
measurements = zeros(ni,2);
baseMeasurement = evalExpr(h, false);
errs = repmat(0.1, ni, 2); % simple error model
for j = 1:ni
    measurements(j,:) = subs(baseMeasurement, t, tvals(j));
end
measurements = measurements + errs; % add error term for measurement

%% Calculate FOCEI obj fun vals
epsi = sym('eps', [ni,2]); % fill these in with dummy values from above
for j = 1:ni
%     hj = subs(h, t, times(j));
    epsi(j,:) = d(j,:) - h;
end
epsi = epsi.'; % nh x nObs matrix

% Objective function value
lij = sym('lij', [ni,1]); % terms in Eq 10
for j = 1:ni
    epsij = subs(epsi(:,j), t, times(j));
    lij(j) = epsij.'*inv(Ri)*epsij + log(det(2*pi*Ri));
end
li = -1/2*sum(lij) - 1/2*eta.'*inv(Omega)*eta - 1/2*log(det(2*pi*Omega));

% Show components
% double(evalExpr(epsi))
% double(evalExpr(Omega))
% double(evalExpr(R))
% double(evalExpr(lij))
% double(evalExpr(-1/2*eta.'*inv(Omega)*eta - 1/2*log(det(2*pi*Omega))))
% double(evalExpr(li))

% Gradient wrt eta
% G = gradient(li, eta);

% Hessian wrt eta - exact
% H = hessian(li, eta);

%% log-likelihood
% logLi = li - 1/2*log(det(-H/(2*pi)));
% objfunval = -2*double(evalExpr(logLi))

%% Gradient of log-likelihood wrt theta
% dlogLF_dtheta = gradient(logLi, theta);
% objfungrad = double(evalExpr(dlogLF_dtheta))

%% Testing individual components in gradient calculation
dOmega_inv_dtheta = sym(zeros(neta,neta,ntheta));
Omega_ = inv(Omega);
for m = 1:ntheta
    for i = 1:neta
        for j = 1:i
            dOmega_inv_dtheta(i,j,m) = diff(Omega_(i,j), theta(m));
            dOmega_inv_dtheta(j,i,m) = diff(Omega_(i,j), theta(m));
        end
    end
end
dOmega_inv_dtheta_vals = double(evalExpr(dOmega_inv_dtheta));

dRi_detai = sym(zeros(neta,neta,neta)); % all zeros in this example
for k = 1:neta
    for ih = 1:nh
        for jh = 1:ih
            dRi_detai(ih,jh,k) = diff(Ri(ih,jh), eta(k));
        end
    end
end
dRi_detai_vals = double(evalExpr(dRi_detai));

dRi_dtheta = sym(zeros(neta,neta,ntheta));
for m = 1:ntheta
    for ih = 1:nh
        for jh = 1:ih
            dRi_dtheta(ih,jh,m) = diff(Ri(ih,jh), theta(m));
            dRi_dtheta(jh,ih,m) = diff(Ri(ih,jh), theta(m));
        end
    end
end
dRi_dtheta_vals = zeros(neta,neta,ntheta,ni);
for j = 1:ni
    dRi_dtheta_vals(:,:,:,j) = double(evalTime(evalExpr(dRi_dtheta), tvals(j)));
end

dBstar_dtheta = sym(zeros(nh,nh,ntheta));
for m = 1:ntheta
    for ih = 1:nh
        for jh = 1:ih
            dBstar_dtheta(:,:,m) = -2*Ri_*dRi_dtheta(:,:,m)*Ri_;
        end
    end
end
dBstar_dtheta_vals = zeros(neta,neta,ntheta,ni);
for j = 1:ni
    dBstar_dtheta_vals(:,:,:,j) = double(evalTime(evalExpr(dBstar_dtheta), tvals(j)));
end

depsi_detai = sym(zeros(nh,neta));
for k = 1:neta
    epsij = epsi(:,j);
    depsi_detai(:,k) = diff(epsij, eta(k));
end
depsi_detai_vals = zeros(nh,neta,ni);
for j = 1:ni
    depsi_detai_vals(:,:,j) = double(evalTime(evalExpr(depsi_detai), tvals(j)));
end

depsi_dtheta = sym(zeros(nh,ntheta));
for m = 1:ntheta
    epsij = epsi(:,j);
    depsi_dtheta(:,m) = diff(epsij, theta(m));
end
depsi_dtheta_vals = zeros(nh,ntheta,ni);
for j = 1:ni
    depsi_dtheta_vals(:,:,j) = double(evalTime(evalExpr(depsi_dtheta), tvals(j)));
end

% Eq 47 can be calculated directly symbolically (is this the exact version of an approximation)
dli_dtheta = gradient(li, theta);
d2li_detaidtheta = sym(zeros(neta,ntheta));
for m = 1:ntheta
    d2li_detaidtheta(:,m) = gradient(dli_dtheta(m), eta);
end
d2li_detaidtheta_vals = double(evalExpr(d2li_detaidtheta));

d2li_detai2 = hessian(li, eta);
d2li_detai2_vals = double(evalExpr(d2li_detai2));

detaistar_dtheta = -inv(d2li_detai2)*d2li_detaidtheta;
detaistar_dtheta_vals = double(evalExpr(detaistar_dtheta));

% Eq 34
depsistar_dtheta = sym(zeros(nh,ntheta));
for m = 1:ntheta
    depsistar_dtheta(:,m) = depsi_dtheta(:,m) + depsi_detai*detaistar_dtheta(:,m);
end
depsistar_dtheta_vals = zeros(nh,ntheta,ni);
for j = 1:ni
    depsistar_dtheta_vals(:,:,j) = double(evalTime(evalExpr(depsistar_dtheta), tvals(j)));
end

% Eq 35 - 2nd term = 0 from problem
dRi_detai = sym(zeros(nh,nh,neta)); % 0 from the problem
dRistar_dtheta = dRi_dtheta;
dRistar_dtheta_vals = dRi_dtheta_vals;

% Eq 37
d_dtheta_depsi_detai

0;

%% Helper functions

    function expr = evalTime(expr, tval)
        expr = subs(expr, t, tval);
    end

    function expr = evalExpr(expr, usesim)
        % Helper function that evaluates/substitutes symbolics -> numbers
        % usesim = false -> don't substitute simulated measurements - for 1st
        %   simulating test data
        if nargin < 2
            usesim = true;
        end
        
        expr = subs(expr, x0, x0_vals);
        expr = subs(expr, k_0, k_0_vals);
        expr = subs(expr, eta, eta_vals); % yhat predict these to be 0 for FO, mode for FOCE
        expr = subs(expr, omega, omega_vals);
        expr = subs(expr, sigma, sigma_vals);
        expr = subs(expr, times, tvals);
        
        if usesim
            expr = subs(expr, d, measurements);
        end
    end

end