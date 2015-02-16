function con = InitialValueExperiment(m, tF, s, inp, dos, name)
%InitialValueExperiment Construct a KroneckerBio experimental conditions
%   structure describing an initial value problem
%
%   con = Experiment(m, tF, s, inp, dos, name)
%
%   Inputs
%   m: [ model struct scalar ]
%       The KroneckerBio model for which these experiments will be run
%   tF: [ nonnegative scalar ]
%       The time at which this experiment ends
%   s: [ nonnegative vector ns {m.s} ]
%       The values of the seed parameters
%   inp: [ input struct scalar | handle @(t) returns nonegative vector nu |
%          nonnegative vector nu {m.u} ]
%       The definition of the input species values
%   dos: [ dose struct scalar {doseZero(m)} ]
%       The definition of the dose amounts and schedule
%   name: [ string {''} ]
%       An arbitrary name for the experiment
%
%   Outputs
%   con: [ experiment struct scalar ]
%       The KroneckerBio experimental conditions structure
%
%   For the meanings of the fields of con see "help Uzero"

% (c) 2015 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

if nargin < 6
    name = [];
    if nargin < 5
        dos = [];
        if nargin < 4
            inp = [];
            if nargin < 3
                s = [];
            end
        end
    end
end

if isempty(s)
    s = m.s;
end
if isempty(inp)
    inp = inputConstant(m, m.u);
end
if isempty(dos)
    dos = doseZero(m);
end
if isempty(name)
    name = '';
end

% m
assert(isscalar(m) && is(m, 'Model'), 'KroneckerBio:Experiment:m', 'm must be a Model')
m = keepfields(m, {'Type', 's', 'u', 'ns', 'nu'});
nu = m.nu;

% tF
assert(isscalar(tF) && isreal(tF) && tF >= 0, 'KroneckerBio:Experiment:tF', 'tF must be a scalar real number greater than or equal to zero')

% s
assert(numel(s) == m.ns, 'KroneckerBio:Experiment:s', 's must a vector with length equal to m.ns')
s = vec(s);

% inp
if isnumeric(inp)
    assert(numel(inp) == m.nu, 'KroneckerBio:Experiment:inp', 'inp, when numeric, must have a length of m.nu')
    inp = inputConstant(m, inp);
end

assert(is(inp, 'Input'), 'KroneckerBio:Experiment:inp', 'inp must be an Input')

% dos
assert(is(dos, 'Dose'), 'KroneckerBio:Experiment:dos', 'dos must be a Dose')

% name
assert(ischar(name), 'KroneckerBio:Experiment:name', 'name must be a string')

% Build experiment
con.Type = 'Experiment:InitialValue';
con.Name = name;
con.tF = tF;
con.nu = m.nu;
con.ns = m.ns;
con.nq = numel(inp.q);
con.nh = numel(dos.h);
con.s  = s;
con.q  = inp.q;
con.h  = dos.h;
con.u  = @(t)inp.u(t,inp.q);
try
    % See if u is vectorized
    test = con.u([1,2]);
    assert(size(test) == [nu,2])
catch
    % Nope, vectorize it
    con.u = @(t)vectorize_u(t, @(ti)inp.u(ti,inp.q));
end
con.dudq = @(t)inp.dudq(t,inp.q);
con.d2udq2 = @(t)inp.d2udq2(t,inp.q);
con.d  = @(t)dos.d(t,dos.h);
con.dddh = @(t)dos.dddh(t,dos.h);
con.d2ddh2 = @(t)dos.d2ddh2(t,dos.h);
con.SteadyState = false;
con.Periodic = false;
con.Discontinuities = vec(unique([inp.discontinuities; dos.discontinuities]));
con.Update = @update;

    function val = vectorize_u(t, f_u)
        nt = numel(t);
        val = zeros(nu,nt);
        for i = 1:nt
            val(:,i) = f_u(t(i));
        end
    end

    function con_out = update(s, q, h)
        con_out = InitialValueExperiment(m, tF, s, inp.Update(q), dos.Update(h), name);
    end
end
