m = equilibrium_model;
% m = equilibrium_model_analytic;


m = m.UpdateField(struct('Name', 'eqm2'));

k = m.k;
knew = k + 2;
m = m.Update(knew);