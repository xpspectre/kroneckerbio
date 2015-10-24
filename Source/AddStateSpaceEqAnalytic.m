function m = AddStateSpaceEqAnalytic(m, ddt, rhs)
% Add a state space equation right hand side to Model.Analytic. Adding these
% expressions to a model is incompatible with adding reactions.
% Inputs:
%   m [ Model.Analytic struct ]
%       The model to which the rule will be added
%   ddt [ string ]
%       Name of state that changes over time
%   rhs [ string ]
%       RHS of system of ODEs, evolution of state in ddt over time
% Inputs:
%   m [ Model.Analytic struct ]
%       The model with the state space RHS added.

% Make sure no reactions are present
assert(isempty(m.Reactions), 'AddStateSpaceEqAnalytic:hasReactions', 'Model contains reactions, which are incompatible with directly specified state space eqs.')

% Add equation as rule
name = ['d' regexprep(ddt, '\W', '_') '_dt']; % is this munging necessary - names can be anything?
m = AddRule(m, name, ddt, rhs, 'rate');
