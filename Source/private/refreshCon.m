function newCon = refreshCon(m, con)
%refreshCon Update Kronecker Bio experimental conditions to a new model
%
%   Whenever either the topology or parameters of a model change, the
%   experimental condition structures should be updated because they may
%   depend on the changes. This function loops over the experimental
%   conditions in order to accomplish this common task.

% (c) 2015 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

nTop = numel(m);
nCon = size(con,1);
newCon = experimentZero([nCon,nTop]);
for iTop = 1:nTop
    for iCon = 1:nCon
        newCon(iCon,iTop) = con(iCon,iTop).Update(con(iCon,iTop).s, con(iCon,iTop).q, con(iCon,iTop).h);
    end
end
