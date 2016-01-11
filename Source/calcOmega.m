function fOmega = calcOmega(Omega, thetas)
% Return Omega, its inverse, its derivative, and its derivative of each
% inverse term as function handles
%
% Inputs:
%   Omega [ neta x neta lower triangular cell matrix of strings ]
%       Cell matrix representing expressions for Omega variances/covariances
%   thetas [ ntheta x 1 cell array of strings ]
%       Cell array containing names of all thetas, including omegas. This
%       is the order in which d_dtheta will be calculated. Note: this also
%       includes all etas here, which will be removed when called in
%       observationFocei.
%
% Outputs:
%   fOmega [ struct ]
%       Struct containing function handles to Omega functions:
%       .Omega [ function handle: (ntheta x 1 double array) -> (neta x neta double matrix) ]
%           Function that returns Omega
%       .OmegaI [ function handle: (ntheta x 1 double array) -> (neta x neta double matrix) ]
%           Function that returns inverse of Omega, inv(Omega)
%       .dOmega_dtheta [ function handle: (ntheta x 1 double array) -> (neta x neta x ntheta double 3-tensor) ]
%           Function that returns derivative of each element of Omega w.r.t.
%           theta, where 3rd dimension is the order of theta given in the input
%       .dOmegaI_dtheta [ function handle: x (ntheta x 1 double array) -> (neta x neta x ntheta double 3-tensor) ]
%           Function that returns derivative of each element of inv(Omega)
%           w.r.t theta, where 3rd dimension is the order of theta given in the input
%
% Note/TODO:
%   - This function uses the symbolic toolbox but returns function handles
%   that shouldn't need it - run this function when building model/objective but
%   only use function handles when fitting
%   - Assumes entries in input omega can be directly turned into symbolics - i.e., they have valid names
%   - Actually enforcing the separation of etas would lead to just as messy
%   of code, so it isn't done for now

% Sanity checks
assert(size(Omega,1) == size(Omega,2), 'calcOmega:NonSquareOmega', 'Omega must be square')

% Clean up theta names
thetas = vec(thetas);
thetas = thetas(~cellfun('isempty',thetas));
assert(iscellstr(thetas), 'calcOmega:InvalidTheta', 'thetas must be a cell array of strings');

ntheta = length(thetas);
th = sym(thetas);

%% Convert Omega into a symbolic symmetric matrix
neta = size(Omega,1);
om = sym(zeros(neta));
for i = 1:neta
    for j = 1:i
        omega_ij = sym(Omega{i,j});
        om(i,j) = omega_ij;
        om(j,i) = omega_ij;
    end
end
omf = matlabFunction(om, 'Vars', {th});

%% Inverse
omI = inv(om);
omIf = matlabFunction(omI, 'Vars', {th});


%% Derivative
dom_dth = sym(zeros(neta,neta,ntheta));
for i = 1:neta
    for j = 1:i
        der_ij = gradient(om(i,j), th);
        dom_dth(i,j,:) = der_ij;
        dom_dth(j,i,:) = der_ij;
    end
end
dom_dthf = matlabFunction(dom_dth, 'Vars', {th});

%% Derivative of inverse
domI_dth = sym(zeros(neta,neta,ntheta));
for i = 1:neta
    for j = 1:i
        der_ij = gradient(omI(i,j), th);
        domI_dth(i,j,:) = der_ij;
        domI_dth(j,i,:) = der_ij;
    end
end
domI_dthf = matlabFunction(domI_dth, 'Vars', {th});


%% Assign final function handles
fOmega = [];
fOmega.Omega          = omf;
fOmega.OmegaI         = omIf;
fOmega.dOmega_dtheta  = dom_dthf;
fOmega.dOmegaI_dtheta = domI_dthf;

end