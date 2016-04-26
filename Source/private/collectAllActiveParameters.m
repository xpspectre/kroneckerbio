function T = collectAllActiveParameters(m, con, ComponentMap, Tlocal2T)
% Assemble parameters to overall Theta vector expected by FitObjective

% Assemble parameters collected from models and conditions into form
%   Tlocal2T expects
nCon = length(con);
[uniqueConInds, ia] = unique(ComponentMap(:,2));
uniqueComponentMap = ComponentMap(ia,:);
assert(nCon == length(uniqueConInds));
Tlocal = cell(nCon,1);
for iCon = 1:nCon
    m_i = m(uniqueComponentMap(iCon,1));
    con_i = con(uniqueComponentMap(iCon,2));
    
    k = m_i.k;
    s = con_i.s;
    q = con_i.q;
    h = con_i.h;
    Tlocal{iCon} = {k,s,q,h}; 
end

T = Tlocal2T(Tlocal);