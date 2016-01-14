function RemoveState(this, name)
% RemoveState Remove a state by name. Name can be an unqualified species name or
% a qualified compartment.species. An unqualified species will be removed if it
% is unique in the model.

import Parser.isFullSpeciesName

if isFullSpeciesName(name)
    existing_full_names = strcat(vec({this.States.Compartment}), '.', vec({this.States.Name}));
    ind = strcmp(name, existing_full_names);
    
    assert(nnz(ind) ~= 0, 'KroneckerBio:Model:RemoveState:StateNotFound', 'State with name %s not found in model', name)
else
    ind = strcmp(name, {this.States.Name});
    
    assert(nnz(ind) ~= 0, 'KroneckerBio:Model:RemoveState:StateNotFound', 'State with name %s not found in model', name)
    assert(nnz(ind) == 1, 'KroneckerBio:Model:RemoveState:AmbiguousStateRemoval', 'Multiple states were found with the name "%s", use the full name ("compartment.%s") to remove it', name, name)
end

this.States(ind,:) = [];
this.nx = this.nx - 1;

this.Ready = false;