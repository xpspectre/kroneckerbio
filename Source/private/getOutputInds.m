function [hInds, RiInds, RiPos] = getOutputInds(yNames, outputs)
% Get indices of output types for NLME fitting. Outputs part of Ri have the Ri__.
% Ri must be a symmetric matrix, whose nonzero lower
% triangular elements must be copied to the top triangle. Only outputs with
% corresponding entries in Ri are part of h - other outputs don't have
% measurements or error models. Sparse - elements missing in Ri that would have
% a corresponding h like a 0 covariance are skipped and assumed 0.
%
% Inputs:
%   yNames [ cell vector of strings ]
%       Names of all outputs in models
%   outputs [ cell vector of strings ]
%       Names of outputs for measurements in objective function in order
%
% Outputs:
%   hInds [ nh x 1 double vector ]
%       Position of ih-th h
%   RiInds [ nh x 1 double vector ]
%       Position of the k-th Ri
%   RiPos [ nh x 2 double vector ]
%       Each row is the (ih,jh)-th position of the corresponding k-th entry
%       in Ri. The ordering of Ri matches the ordering of h

% Get indices corresponding to h and Ri
ny = length(yNames);
RiInds  = zeros(ny,1);
hInds    = zeros(ny,1);
for i = 1:ny
    yName = yNames{i};
    if ismember(yName, outputs);
        hInds(i) = i;
    elseif ~isempty(regexp(yName, '^Ri__', 'once'))
        RiInds(i) = i;
    end % ignore extraneous non-h, non-Ri outputs
end
RiInds(RiInds==0) = [];
hInds(hInds==0)     = [];

% Get (ih,jh) positions of Ri
nRi = length(RiInds);
RiPos = zeros(nRi,2);
RiNames = yNames(RiInds);
for i = 1:nRi
    tokens = regexp(RiNames{i}, '^Ri__(.+)__(.+)$', 'tokens', 'once')';
    assert(iscell(tokens) && length(tokens) == 2, 'observationNLME:getOutputInds:invalidRiTokens', 'Parser didn''t find correct params Ri name %s is referring to', RiNames{i})
    ix = find(ismember(outputs, tokens{1}));
    jx = find(ismember(outputs, tokens{2}));
    RiPos(i,:) = [ix,jx];
end
end