function Finalize(this, opts)
%FinalizeModel Update the mathematical components of the model to reflect
%   changes made to the model.
%
%   This function contains generic model processing steps. Some steps come from
%   the common Model class while others come from the appropriate subclasses.

if nargin < 2
    opts = [];
end

this.qualifyNames;
this.extractSpecies;
this.calculateDerivatives(opts);

