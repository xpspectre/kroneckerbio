function [thetaMask, etaMask] = getParameterMask(names)
% Return logical matrices corresponding to set of thetas and set of
%   etas for outer and inner optimization, respectively, based on
%   parameter name prefixes.
%
% Inputs:
%   names [ cell array of strings ]
%       Parameter names
%
% Outputs:
%   thetaMask [ length(names) x 1 logical vector ]
%       Positions of thetas
%   etaMask [ length(names) x 1 logical vector ]
%       Positions of etas

etaFound = regexp(names, '^eta__', 'once');
thetaMask = cellfun(@isempty, etaFound);
etaMask = ~thetaMask;
