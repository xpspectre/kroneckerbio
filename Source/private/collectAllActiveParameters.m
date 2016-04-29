function [T, details] = collectAllActiveParameters(m, con, ComponentMap, Tlocal2T)
% Assemble parameters to overall Theta vector expected by FitObjective
%
% Inputs:
%	m [ nModel x 1 model struct vector ]
%       The models that will be simulated
% 	con [ nCon x 1 experiment struct vector ]
%       The experimental conditions under which the model will be simulated
%   ComponentMap [ nObjx 3 double matrix ]
%       Model+condition+objective mapping with 1 row for each objective from fit
%       options
%   Tlocal2T [ function handle (nCon x 1 cell vector of 4 x 1 cell vectors
%           of nXi x 1 double vectors) -> (nT x 1 double vector) ]
%      Function from ParamMapper that converts format pulled from
%      model+condition to overall Theta vector.
%
% Outputs:
%   T [ nT x 1 double vector ]
%       Overall Theta vector of parameteres to fit from FitObjective
%    details [ nT x [Name,Type,Ind] Table ]
%        Table of parameter name, type {'k','s','q','h'} and the
%        model (for k) or experiment (for s,q,h) it came from
%

% Assemble parameters collected from models and conditions into form
%   Tlocal2T expects
nCon = length(con);
Tlocal = cell(nCon,1);
for iCon = 1:nCon
    m_i = m(ComponentMap{iCon,1});
    con_i = con(ComponentMap{iCon,2});
    
    k = m_i.k;
    s = con_i.s;
    q = con_i.q;
    h = con_i.h;
    Tlocal{iCon} = {k,s,q,h}; 
end

[T, details_] = Tlocal2T(Tlocal);

% Process details_ into nicer details if requested
if nargout == 2
    paramTypes = {'k','s','q','h'};
    nT = numel(T);
    
    Name = cell(nT,1);
    Type = cell(nT,1);
    Ind = cell(nT,1);
    
    for i = 1:nT
        Ind{i} = details_(i,3);
        Type{i} = paramTypes{details_(i,1)};
        switch Type{i}
            case 'k'
                kNames = {m(Ind{i}).Parameters.Name};
                Name{i} = kNames{details_(i,2)};
            case 's'
                parentModelInd = [ComponentMap{:,2}] == Ind{i};
                sNames = {m(parentModelInd).Seeds.Name};
                Name{i} = sNames{details_(i,2)};
            otherwise
                error('Names for q and h not implemented yet')
        end
    end
    
    details = table(Name,Type,Ind);
end