function Tlocal = mapT2Tlocal(T, Tmap)
% Map combined T vector of fit params for 1 condition to standard form.
% Input:
%   T [ nT x 1 double vector ]
%       Format of combined fit params for this condition
%   Tmap [ 1 x 4 cell vector of nX x 1 logical vectors ]
%       Format of which params are used
% Output:
%   Tlocal [ 1 x 4 cell vector of nX x 1 double vectors ]
%       Format of this condition's fit params
Tlocal = cell(1,4);
Tidx = 0;
for i = 1:4
    Tstart = Tidx + 1;
    Tend = Tstart + nnz(Tmap{i}) - 1;
    Tsub = T(Tstart:Tend);
    
    Tlocal_i = zeros(length(Tmap{i}),1);
    Tlocal_i(Tmap{i}) = Tsub;
    Tlocal{i} = Tlocal_i;
    
    Tidx = Tend;
end