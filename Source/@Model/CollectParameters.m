function [T, info] = CollectParameters(this)
% Return parameters of model in order they appear
%
% Outputs:
%   T [ this.nk x 1 vector of doubles ]
%       Values of the parameters
%   info [ this.nk x 1 struct array ]
%       Detailed info on parameters, with fields:
%       .Type [ string ]
%           Parameter type, kinetic parameters 'k' in models for now
%       .Name [ string ]
%           Parameter name

T = vec([this.Parameters.Value]);

% Also return info
if nargout > 1
    info = struct('Type', 'k' ,'Name', {this.Parameters.Name});
end
