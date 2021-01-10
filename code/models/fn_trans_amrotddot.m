function output = fn_trans_amrotddot(coeffs,eqdata)
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

%% RJW Lift (includes deviation) - NOT FITTING AM COEFFS HERE, BUT WOULD BE
%% WORTH DOING SO
% Unpack Coefficients:
Clmax = coeffs(1);
Camrotddot = coeffs(2);

% Compute Force
wx = defdot-strokedot.*sin(dev);
wxdot = defddot - (strokeddot.*sin(dev)+strokedot.*devdot.*cos(dev));
wz = devdot.*cos(def) + strokedot.*sin(def).*cos(dev); % Vo = r_i*omega_z
wy = -strokedot.*cos(def).*cos(dev) + devdot.*sin(def); % Wo = r_i*omega_y
wydot = -strokeddot.*cos(def).*cos(dev)+strokedot.*defdot.*sin(def).*cos(dev)+...
    strokedot.*devdot.*cos(def).*sin(dev)+devddot.*sin(def)+devdot.*defdot.*cos(def);
wh = sqrt(wz.^2+wy.^2); % angular velocity of the wing hinge
alpha = atan2(-wy,wz); % angle of attack relative to instantaneous velocity
beta = alpha+def; % angle between velocity vector and vertical reference
Cl = Clmax.*sin(2*alpha); % coefficient of lift
Lrjw = .5*rho.*chord*(1/3).*(outerspan.^3-innerspan.^3).*Cl.*wh.^2.*sin(beta);
% lambdaz = pi*rho*(chord/2).^2; % kg/m^3 * m^2
lambdazw = -pi*rho*(chord/2).^2.*(chord/2-upperchord); % kg/m^3 * m^2 * m
% Z01 = -lambdaz*.5.*(wydot-wx.*wz).*-(outerspan.^2-innerspan.^2);
Z02 = Camrotddot*-lambdazw.*wxdot.*(outerspan-innerspan);
% Z0spanwise = Z01+Z02;
% AMrjw = Z0spanwise.*-sin(def);
% output = Lrjw+AMrjw+inertial+inerttheta;
% AMTrans = Z01.*-sin(def);
AMRotdot = Z02.*-sin(def);
output = Lrjw+AMRotdot+inertial+inerttheta;