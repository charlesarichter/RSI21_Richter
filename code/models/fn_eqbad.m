function output = fn_eqbad(coeffs,eqdata)
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
C4 = coeffs(4);
C5 = coeffs(5);
C6 = coeffs(6);
C7 = coeffs(7);
C8 = coeffs(8);
C9 = coeffs(9);

output = (14.9488 + 13.5087*def.*def - 1.38486.*def)./...
  (407.258 + devddot + 17.8828.*devdot +...
  cos(2.69241 - 42.9085./(13.5087.*def.*def) - 163.587.*def)) - 0.0381286;

% Coeffs = [14.9488, 13.5087, 1.38486, 407.258, 17.8828, 2.69241, 42.9085, 13.5087, 163.587, 0.0381286]