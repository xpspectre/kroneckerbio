function [thetaInds, etaInds] = getParameterInds(names)
% Get indices of parameter types for NLME fitting. omegaInds are a subset of thetaInds and
%   included separately for convenience. Omega has the same order as Eta.
%   Any non-Eta parameter is in Omega.
%
% Inputs:
%   names [ cell array of strings ]
%       Parameter names
%
% Outputs:
%   thetaInds [ ntheta x 1 double vector ]
%       Position of the m-th theta
%   etaInds [ neta x 1 double vector ]
%       Position of the k-th eta

[thetaMask, etaMask] = getParameterMask(names);
thetaInds = find(thetaMask);
etaInds = find(etaMask);
