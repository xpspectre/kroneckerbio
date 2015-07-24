function X = matSqrtInv(C)
% Calculates the matrix inverse square root
% Mathematically:
%   Given a nonsingular matrix C in Cnxn, a matrix X s.t. X^2*C = I is called
%   the matrix inverse square root.
% Implementation from Min-Gul Kim, et al., 2015; further description in Sherif,
% 1991
[V,D] = eig(C);
d = abs(diag(D));
d2 = diag(1./sqrt(d));
X = V*d2*V';
end