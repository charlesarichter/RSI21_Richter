function [output,term1,term2,term3] = fn_eq4(coeffs,eqdata)
%% Unpack Data
force = eqdata(:,1);
stroke = eqdata(:,2);strokedot = eqdata(:,3);strokeddot = eqdata(:,4);
def = eqdata(:,5);defdot = eqdata(:,6);defddot = eqdata(:,7);
dev = eqdata(:,11);devdot = eqdata(:,12);devddot = eqdata(:,13);
% dev = zeros(size(dev));devdot = zeros(size(devdot));devddot = zeros(size(devddot));
% span = 1e-3*eqdata(:,17);innerspan = 1e-3*eqdata(:,18);outerspan = 1e-3*eqdata(:,19);
% chord = 1e-3*eqdata(:,22);lowerchord = 1e-3*eqdata(:,20);upperchord = 1e-3*eqdata(:,21);
span = eqdata(:,17);innerspan = eqdata(:,18);outerspan = eqdata(:,19);
chord = eqdata(:,22);lowerchord = eqdata(:,20);upperchord = eqdata(:,21);
Xcm	= eqdata(:,23);Ycm	= eqdata(:,24);Zcm= eqdata(:,25);
Ixx	= eqdata(:,26);Ixy	= eqdata(:,27);Ixz= eqdata(:,28);
Iyx	= eqdata(:,29);Iyy	= eqdata(:,30);Iyz= eqdata(:,31);
Izx	= eqdata(:,32);Izy	= eqdata(:,33);Izz= eqdata(:,34);

%% Compute Output
C1 = coeffs(1);
C2 = coeffs(2);
C3 = coeffs(3);

term1 = C1.*span.*Ycm.*defdot.*defdot.*cos(def);
term2 = C2.*strokedot.*span.*Xcm.*def;
term3 = -C3.*span.*devddot;
output = term1 + term2 + term3;

% output = C1.*span.*Ycm.*defdot.*defdot.*cos(def)...
%     + C2.*strokedot.*span.*Xcm.*def...
%     - C3.*span.*devddot;

% Coeffs = [1.4526972e-8 1.2489714e-7 1.0872654e-6]
