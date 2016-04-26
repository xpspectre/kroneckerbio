function [m, con] = updateAllActiveParameters(m, con, T, ComponentMap, T2Tlocal)
% Update individual models and conditions with overall Theta vector from
%   FitObjective
%
% Inputs:
%	m [ nModel x 1 model struct vector ]
%       The models that will be simulated
% 	con [ nCon x 1 experiment struct vector ]
%       The experimental conditions under which the model will be simulated
%   T [ nT x 1 double vector ]
%       Overall Theta vector of parameteres to fit from FitObjective
%   ComponentMap [ nObjx 3 double matrix ]
%       Model+condition+objective mapping with 1 row for each objective from fit
%       options
%   T2Tlocal [ function handle (nT x 1 double vector) -> (nCon x 1 cell
%           vector of 4 x 1 cell vectors of nXi x 1 double vectors) ]
%       Function from ParamMapper that converts overall Theta vector to format
%       each model+condition can use
%
% Outputs:
%	m [ nModel x 1 model struct vector ]
%       The models with updated k
% 	con [ nCon x 1 experiment struct vector ]
%       The experimental conditions with updated s,q,h

nCon = length(con);
[uniqueConInds, ia] = unique(ComponentMap(:,2));
uniqueComponentMap = ComponentMap(ia,:);
assert(nCon == length(uniqueConInds));

Tlocals = T2Tlocal(T);
for iCon = 1:nCon
    m_i = m(uniqueComponentMap(iCon,1));
    con_i = con(uniqueComponentMap(iCon,2));
    
    % Combine with existing params existing params
    Tlocalold = {m_i.k, con_i.s, con_i.q, con_i.h};
    Tlocal = combineCellMats(Tlocalold, Tlocals{iCon});
    
    [k,s,q,h] = deal(Tlocal{:});
    
    m(uniqueComponentMap(iCon,1)) = m_i.Update(k);
    con(uniqueComponentMap(iCon,2)) = con_i.Update(s,q,h);
end
