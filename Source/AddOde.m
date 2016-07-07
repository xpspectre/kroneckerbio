function m = AddOde(m, ddt, rhs)
% Add an ODE right hand side directly to an analytic model. Adding these
%   expressions to a model is incompatible with adding reactions.
%
% Inputs:
%   m [ Model.Analytic struct ]
%       The model to which the rule will be added
%   ddt [ string ]
%       Name of state that changes over time
%   rhs [ string ]
%       RHS of system of ODEs, evolution of state in ddt over time
%
% Outputs:
%   m [ Model.Analytic struct ]
%       The model with the state space RHS added.
%
% Note: This functionality is implemented using a new type of rule.

% Make sure no reactions are present
assert(isempty(m.Reactions), 'AddOde:HasReactions', 'Model contains reactions, which are incompatible with directly specified ODEs.')

% Add equation as rule
name = ['d(' ddt ')/dt'];
m = AddRule(m, name, ddt, rhs, 'rate');
