function outputs = fixOutputs(outputs)
% Convert output argument to right form for NLME observations
%
% Inputs:
%   outputs [ string | nOutputs cell array of strings ]
%
% Outputs:
%   outputs [ nOutputs x 1 cell vector of strings ]

if iscellstr(outputs)
    outputs = vec(outputs);
elseif ischar(outputs)
    outputs = {outputs};
else
    error('KroneckerBio:fixOutputs:InvalidOutputArg', 'Output must be a string or cell array of strings')
end
