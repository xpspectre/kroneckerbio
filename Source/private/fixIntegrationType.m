function [continuous, complex, tGet] = fixIntegrationType(con, obj)
%FIXINTEGRATIONTYPE Standardize continuous, complex, and tGet
%
%   There are two attributes of objective functions that affect how the
%   system is integrated. Continuous indicates that there is a continuous
%   aspect to the objective function. If any objective function has this
%   atribute, then the integration must integrate it. Complex indicates
%   that the objective must have the complete solution in order to compute
%   the objective (like to optimize the peak time or peak value of a
%   species). Objective functions that are not complex must have
%   discreteTimes indicated to tell when the objective function needs to be
%   evaluated.

% (c) 2010 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

% Constants
[nCon, nObj] = size(obj);

% Initialize
continuous = false(nCon,1);
complex = false(nCon,1);
tGet = cell(nCon,1);

% Loop through each objective to determine if any require special treatment
for iCon = 1:nCon
    for iObj = 1:nObj
        continuous(iCon) = continuous(iCon) || obj(iCon,iObj).Continuous;
        complex(iCon) = complex(iCon) || obj(iCon,iObj).Complex;
    end
    
    % Only for non-complex integrations is tGet useful
    if ~complex(iCon)
        discreteTimes = [];
        for iObj = 1:nObj
            discreteTimes = [discreteTimes; vec(obj(iCon,iObj).DiscreteTimes)];
        end
        if continuous(iCon)
            % We need the final point as well
            discreteTimes = [discreteTimes; con(iCon).tF];
        end
        tGet{iCon} = unique(discreteTimes);
    end
end