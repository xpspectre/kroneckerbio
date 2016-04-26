function [m, con] = updateAllActiveParameters(m, con, T, ComponentMap, T2Tlocal)
% Update individual models and conditions with overall Theta vector from
%   FitObjective

nCon = length(con);
[uniqueConInds, ia] = unique(ComponentMap(:,2));
uniqueComponentMap = ComponentMap(ia,:);
assert(nCon == length(uniqueConInds));
Tlocal = T2Tlocal(T);
for iCon = 1:nCon
    m_i = m(uniqueComponentMap(iCon,1));
    con_i = con(uniqueComponentMap(iCon,2));
    
    k = Tlocal{iCon}{1};
    s = Tlocal{iCon}{2};
    q = Tlocal{iCon}{3};
    h = Tlocal{iCon}{4};
    
    m(uniqueComponentMap(iCon,1)) = m_i.Update(k);
    con(uniqueComponentMap(iCon,2)) = con_i.Update(s,q,h);
end
1;