function output = fn_planarlift_planarinertial(coeffs,eqdata)
%% Unpack Data
force = eqdata(:,1);
stroke = eqdata(:,2);strokedot = eqdata(:,3);strokeddot = eqdata(:,4);
def = eqdata(:,5);defdot = eqdata(:,6);defddot = eqdata(:,7);
dev = eqdata(:,11);devdot = eqdata(:,12);devddot = eqdata(:,13);
% dev = zeros(size(dev));devdot = zeros(size(devdot));devddot = zeros(size(devddot));
span = 1e-3*eqdata(:,17);innerspan = 1e-3*eqdata(:,18);outerspan = 1e-3*eqdata(:,19);
chord = 1e-3*eqdata(:,22);lowerchord = 1e-3*eqdata(:,20);upperchord = 1e-3*eqdata(:,21);
Xcm	= eqdata(:,23);Ycm	= eqdata(:,24);Zcm= eqdata(:,25);
Ixx	= eqdata(:,26);Ixy	= eqdata(:,27);Ixz= eqdata(:,28);
Iyx	= eqdata(:,29);Iyy	= eqdata(:,30);Iyz= eqdata(:,31);
Izx	= eqdata(:,32);Izy	= eqdata(:,33);Izz= eqdata(:,34);

%% Other Constants
rho720 = 1.2159e+003; % 1.092*(1/1000)*(100/1)^3; % kg/m^3
thickness = 9.6000e-005; % m
rho = 1.2; % kg/m^3, density of air

%% Inertial
inertial = -(chord/2-upperchord).*chord*thickness.*span*rho720.*(defddot.*sin(def) + defdot.^2.*cos(def));
inerttheta = -.5*(innerspan+outerspan).*span.*chord*thickness*rho720.*(devddot.*cos(dev) - devdot.^2.*sin(dev));

%% Planar U-Squared Lift, Planar Inertial Forces - COMPLEXITY 47
% Unpack Coefficients:
Clmax = coeffs(1);

% Compute Force
beta = atan2(strokedot,0);
alpha = beta-def;
Cl = Clmax*sin(2*alpha);
lift = .5*rho.*chord*(1/3).*(outerspan.^3-innerspan.^3).*Cl.*strokedot.^2.*sin(beta);
output = lift+inertial;