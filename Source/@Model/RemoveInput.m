function RemoveInput(this, name)
% RemoveInput Remove an input by name. Name can be an unqualified species name or
% a qualified compartment.species. An unqualified species will be removed if it
% is unique in the model.

import Parser.isFullSpeciesName

if isFullSpeciesName(name)
    existing_full_names = strcat(vec({this.Inputs.Compartment}), '.', vec({this.Inputs.Name}));
    ind = strcmp(name, existing_full_names);

    assert(nnz(ind) ~= 0, 'KroneckerBio:Model:RemoveInput:InputNotFound', 'Input with name %s not found in model', name)
else
    ind = strcmp(name, {this.Inputs.Name});
    
    assert(nnz(ind) ~= 0, 'KroneckerBio:Model:RemoveInput:InputNotFound', 'Input with name %s not found in model', name)
    assert(nnz(ind) == 1, 'KroneckerBio:Model:RemoveInput:AmbiguousInputRemoval', 'Multiple inputs were found with the name "%s", use the full name ("compartment.%s") to remove it', name, name)
end

this.Inputs(ind,:) = [];
this.nu = this.nu - 1;

this.Ready = false;