function RemoveReaction(this, name)
% RemoveReaction Remove a reaction by name. Since reactions don't have to have
% names, this obviously only works when you give them useful names.

% Find added instances of this name
ind = strcmp(name, {this.Reactions.Name});
assert(sum(ind)==1, 'KroneckerBio:Model:RemoveReaction:ReactionNotFound', 'Reaction with name %s not found in model', name)

% Remove all mention of this seed
this.Reactions(ind,:) = [];
this.nr = this.nr - 1;

this.Ready = false;
