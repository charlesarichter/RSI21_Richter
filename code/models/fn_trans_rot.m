function output = fn_trans_rot(coeffs,eqdata)
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

%% Trans_Rot: Complexity ZJW - 12
Ct = coeffs(1); % 1.2
Cr = coeffs(2); % pi; 1.9254;
% m22 = pi*rho*(chord/2).^2;

% wx = defdot-strokedot.*sin(dev);
wz = devdot.*cos(def) + strokedot.*sin(def).*cos(dev);
wy = -strokedot.*cos(def).*cos(dev) + devdot.*sin(def);
% wydot = -strokeddot.*cos(def).*cos(dev)+strokedot.*defdot.*sin(def).*cos(dev)+...
%     strokedot.*devdot.*cos(def).*sin(dev)+devddot.*sin(def)+devdot.*defdot.*cos(def);
a = chord/2;
R2 = (outerspan.^2-innerspan.^2)./2;
R3 = (outerspan.^3-innerspan.^3)./3;

Fxprime = -rho*-wy.*(-2*Ct.*a.*R3.*wz.*-wy./sqrt(wz.^2+wy.^2)+2*Cr*a.^2.*R2.*defdot);
Fyprime = rho*wz.*(-2*Ct.*a.*R3.*wz.*-wy./sqrt(wz.^2+wy.^2)+2*Cr*a.^2.*R2.*defdot);

Lzjw = Fxprime.*cos(def)-Fyprime.*sin(def);
output = Lzjw+inertial+inerttheta;