% Analytic solution and explicit NLME expressions for testing/comparison for
% Theophylline model
% From Min-Gul Kim, et al., 2015
% Note that positions of k = kel and V are switched
function theo_analytic
theoph = load_theo();

% Same starting params as test case
Theta = [1.5 0.15 25];
Eta = [0 0 0];
OM = [0.2 0.1 0.1; 0.1 0.2 0.1; 0.1 0.1 0.2];
SG = [0.1 0; 0 0.1];

nPatients = height(theoph);
vals = zeros(nPatients,1);
for i = 1:nPatients
    data = theoph(i,:);
    [Fi, Gi, Hi] = pred(Theta, Eta, data);
    Yi = data.Conc{1};
    
    RES = Yi - Fi;
    COV = Gi*OM*Gi' + diag(diag(Hi*SG*Hi'));
    WRES = matSqrtInv(COV)*RES;
    
    vals(i) = log(det(COV)) + WRES'*WRES;
end
val = sum(vals);

end

function [Fi, Gi, Hi] = pred(Theta, Eta, Data)
% Per-patient prediction
dose = 320;
time = Data.Time{1};

ka  = Theta(1)*exp(Eta(1));
kel = Theta(2)*exp(Eta(2));
V   = Theta(3)*exp(Eta(3));

Fi = dose/V * ka/(ka - kel) * (exp(-kel*time) - exp(-ka*time));

G1 = ((dose/V/(ka - kel) - dose/V .* ka/(ka - kel).^2) .* (exp(-kel*time) - exp(-ka*time)) + ...
    dose/V .* ka/(ka - kel) .* (exp(-ka*time).*time)) .* ka;

G2 = (dose/V .* ka/(ka - kel)^2 .* (exp(-kel*time) - exp(-ka*time)) - ...
    dose/V .* ka/(ka - kel) .* (exp(-kel*time).*time)) .* kel; % G3 in paper

G3 = -(dose/V^2 .* ka/(ka - kel) .* (exp(-kel*time) - exp(-ka*time))) .* V; % G2 in paper

H1 = Fi;

H2 = ones(size(Fi));

Gi = [G1, G2, G3];
Hi = [H1, H2];

end