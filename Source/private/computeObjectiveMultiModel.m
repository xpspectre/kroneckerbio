function [G, D] = computeObjectiveMultiModel(m, con, obj, opt)

nTop = length(m);
nT = length(opt.LowerBound); % arbitrary but contains the right number

% Make array of opts, each referring to the corresponding model
opts = fixOpts(m, opt);

if nargout == 1
    G = 0;
    for i = 1:nTop
        G = G + computeObj(m(i), con(i), obj(i,i), opts(i));
    end
end

if nargout == 2
    G = 0;
    D = zeros(nT,1);
    for i = 1:nTop % parfor, but may be overhead in transferinig obj funs
        [Gi, Dii] = computeObjGrad(m(i), con(i), obj(i,i), opts(i));
        
        % Reshape D to match overall T format, placing params in right positions
        %   When a param is shared, all Dii sum their contributions
        Di = mapT2Ti(Dii, i, opt);
        
        G = G + Gi;
        D = D + Di;
    end
end
