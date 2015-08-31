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
%
% These terms are implemented as a special type of rule, which is handled
% specially by finalizeModelAnalytic.

% Increment counter
nz = m.add.nz + 1;
m.add.nz = nz;
m.add.Rules = growRulesAnalytic(m.add.Rules, m.add.nz);

% Make sure no reactions are present
if ~isempty(m.Reactions) || ~isempty(m.add.Reactions)
    error('AddStateSpaceEqAnalytic:hasReactions', 'Model contains reactions, which are incompatible with directly specified state space eqs.')
end

% Make name
name = ['d' regexprep(ddt, '\W', '_') '_dt'];

% Add item
m.add.Rules(nz).Name       = name;
m.add.Rules(nz).ID         = genUID;
m.add.Rules(nz).Target     = ddt;
m.add.Rules(nz).Expression = fixRuleExpression(rhs);
m.add.Rules(nz).Type       = 'state space eq';

m.Ready = false;