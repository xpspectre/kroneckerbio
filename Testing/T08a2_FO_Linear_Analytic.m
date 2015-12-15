% Analytic/symbolic version of simple linear model

function T08a2_FO_Linear_Analytic
clear; close all; clc
rng('default');

% Model

x0 = sym('x0');
k0 = sym('k0');
eta_k = sym('eta_k');
t = sym('t');
err = sym('err');

h = x0 + k0*exp(eta_k)*t;
y = h + err;

omega__k__k = sym('omega__k__k');
Omega = omega__k__k;

sigma__y = sym('sigma__y');
R = sigma__y.^2;

eta = eta_k;
neta = length(eta);
theta = [k0; omega__k__k; sigma__y];
ntheta = length(theta);

% Experimental observations
nObs = 5;
t0 = 0;
tf = 10;
tvals = linspace(0,tf,nObs)';
dobs = repmat(y, nObs, 1);
dobs = subs(dobs, x0, 1);
dobs = subs(dobs, k0, 1.5);
dobs = subs(dobs, eta, 0);

% errs = normrnd(0, sqrt(0.05), nObs, 1);
errs = repmat(0.1, nObs, 1);

% times = sym(linspace(0,tf,nObs));
times = sym('t', [nObs,1]);
epsi = sym(zeros(nObs,1));
% d = sym(zeros(nObs,1));
d = sym('d', [nObs,1]);
for j = 1:nObs
    yj = subs(y, t, times(j));
    epsi(j) = d(j) - yj;
    
    dobs(j) = subs(dobs(j), t, tvals(j));
    dobs(j) = subs(dobs(j), err, errs(j));
end
double(dobs);

% Objective function value
lij = sym('lij', [nObs,1]); % terms in Eq 10
for j = 1:nObs
    lij(j) = epsi(j).'*inv(R)*epsi(j) + log(det(2*pi*R));
end
li = -1/2*sum(lij) - 1/2*eta.'*inv(Omega)*eta - 1/2*log(det(2*pi*Omega));
double(evalExpr(epsi));
double(evalExpr(Omega));
double(evalExpr(R));
double(evalExpr(lij));
double(evalExpr(-1/2*eta.'*inv(Omega)*eta - 1/2*log(det(2*pi*Omega))));
double(evalExpr(li));

% Gradient wrt eta
G = gradient(li, eta);

% Hessian wrt eta - exact
H = hessian(li, eta);

%% Hessian wrt eta - Almquist eq 14
% deps_deta = sym('deps_deta', [1,neta]);
% a = sym('a', [1,neta]);

%% log-likelihood
logLi = li - 1/2*log(det(-H/(2*pi)));
-2*double(evalExpr(logLi))

%% Gradient of log-likelihood wrt theta
% params 1,3,4 (2 is eta) in full system
dlogLF_dtheta = gradient(logLi, theta);
double(evalExpr(dlogLF_dtheta))

0;

%% Helper functions

    function expr = evalExpr(expr)
        % Helper function that evaluates/substitutes symbolics -> numbers
        expr = subs(expr, x0, 1);
        expr = subs(expr, k0, 1.5);
        expr = subs(expr, times, tvals);
        expr = subs(expr, eta, 0); % yhat predict these to be 0 for FO, mode for FOCE
        expr = subs(expr, omega__k__k, 0.2);
        expr = subs(expr, sigma__y, sqrt(0.05));
        expr = subs(expr, d, dobs);
        expr = subs(expr, err, 0); % yhat predicts these to be 0
    end

end