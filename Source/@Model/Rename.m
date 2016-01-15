function Rename(this, name)
% Rename the model to name
%
% Inputs:
%   name [ string ]
%       New model name
%
% Note: I found a call to this in LoadModelMassAction but it wasn't
% previously defined anywhere. I'm adding it because it's useful and easy
% (easier than before, where the closure needs to be updated).

assert(ischar(name), 'KroneckerBio:Model:Rename:InvalidName', 'name must be a string.')
this.Name = name;