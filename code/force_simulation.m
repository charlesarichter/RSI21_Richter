%% Load Data

clear,clc,close all
%----
##folderName = '/Users/charlesrichter/Desktop/CR Shared Drive/';
##spreadsheet = dlmread([folderName 'selecteddata_defcorr_all_70'],',',2,0);
folderName = '/home/crichter/ccsl/';
spreadsheet = dlmread('/home/crichter/ccsl/selecteddata_defcorr_all_70.csv',',',2,0);
% ----
%22,32/33,38,42,51,56,58,61,85,105
% withheld = [22,32,38,42,51,56,58,61,85,105];
% 32,51,58

% Capture only the withheld experiments for testing
exps = [8 26 33 37 39 48 52 53 95 103]; % 1,78 Removed
EXPERIMENT_USED_TO_SHOWCASE_EQ4_COMPONENTS = 16; %%%%%%%%%%%


% totalexps = unique(spreadsheet(:,1));
% exps = totalexps(totalexps>0);

% nExps = length(exps);
% eqdataExps = cell(nExps,1);
% whichExp = cell(nExps,1);
% forceExps = cell(nExps,1);
% for i = 1:nExps
%     eqdataExps{i} = spreadsheet(spreadsheet(:,1)==exps(i),2:end);
%     forceExps{i} = eqdataExps{i}(:,1);
% end
% eqdata = eqdataExps{1};

% 19, 46 - use 41 instead of 46

totalexps = unique(spreadsheet(:,1));
exps = totalexps(totalexps>0);
nExps = max(exps);
eqdataExps = cell(nExps,1);
for i = 1:length(exps)
    whichExp = unique(spreadsheet(spreadsheet(:,1)==exps(i),1));
    eqdataExps{whichExp} = spreadsheet(spreadsheet(:,1)==whichExp,2:end);
end


% figure(4)
% for bla2 = 1:length(exps)
%     disp(exps(bla2))
%     plot(eqdataExps{exps(bla2)}(:,1));
%     pause
% end


%%
% use 52 53 33 48
% use 33 39 53 103? 37?
% toplot = [95 53 33 48 16]; % for exp 94 stuff
% toplot = [39 48 52 53];
% toplot = [8 26 33 37];
toplot = [16 33 48 39 16];
%16 39 48 33

% For exp 95 stuff
% 37 48 52 53
% toplot = [95 37 48 52 16];
% toplot = [22,38,56,58,61,85,105];

% For exp 84 stuff: 33 39 48
% toplot = [84 33 39 48 16];
% toplot = [94 33 37 52];
% [8 26 33 37 39 48 52 53 95 103]; % 1,78 Removed
for bla = 1:length(toplot)
eqdata = eqdataExps{toplot(bla)};

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

%% Plotting Constants
domain = 0:length(force)-1;
axbuf = 1.2;
shadecolor = .9*[1 1 1];

%% Constants
% leadingedge = .010;trailingedge = -.050; % m
% chord = leadingedge-trailingedge; % m
% innerspan = .050;outerspan = .150; % m
% span = outerspan-innerspan; % m
% dspan = .001;dchord = .001; % m
% r = innerspan+dspan/2:dspan:outerspan-dspan/2; % m, Radial blade elements
% ch = trailingedge+dchord/2:dchord:leadingedge-dchord/2; % m, Radial chord elements sum(L,2);
% nBlades = length(r);nChord = length(ch);

thickness = .000096; % m
rho720 = 1215.92; % kg/m^3
rho = 1.2; % kg/m^3, density of air
a = chord/2; % m
yh = a-upperchord; % m
xh = (innerspan+outerspan)/2; % m
lambdaz = pi*rho*a.^2; % kg/m^3 * m^2
lambdazw = -pi*rho*a.^2.*yh; % kg/m^3 * m^2 * m

%% Inertial Reaction
inertial = -yh.*chord.*thickness.*span.*rho720.*(defddot.*sin(def) + defdot.^2.*cos(def));
inertialdev = -xh.*span.*chord.*thickness.*rho720.*(devddot.*cos(dev) - devdot.^2.*sin(dev));

% Also calculated with full integration written out:
inertial2 = .5*(upperchord.^2-lowerchord.^2).*thickness.*span.*rho720.*(defddot.*sin(def) + defdot.^2.*cos(def));
inertialdev2 = -.5*(outerspan.^2-innerspan.^2).*chord.*thickness.*rho720.*(devddot.*cos(dev) - devdot.^2.*sin(dev));

%----
% figure(5),hold on
% ylimit = axbuf*[min(force) max(force)];
% cross = find([0;diff(sign(strokedot))]);
% fill([cross(1) cross(2)-1 cross(2)-1 cross(1)],100*[-1 -1 1 1],shadecolor,'edgecolor',shadecolor)
% fill([cross(3) cross(4)-1 cross(4)-1 cross(3)],100*[-1 -1 1 1],shadecolor,'edgecolor',shadecolor)
% set(gcf,'defaultlinelinewidth',1)
%----

% plot(force)
% plot(-yh.*chord.*thickness.*span.*rho720.*(defdot.^2.*cos(def)),'b')
% plot(-yh.*chord.*thickness.*span.*rho720.*(defddot.*sin(def)),'r')
% plot(-xh.*span.*chord.*thickness.*rho720.*(devddot.*cos(dev)),'g')
% plot(-xh.*span.*chord.*thickness.*rho720.*(-devdot.^2.*sin(dev)),'c')
% 
% plot(1.419276e-8.*1000.*span.*Ycm.*defdot.*defdot.*cos(def),'b.')
% plot(1.2280232e-7*1000.*strokedot.*span.*Xcm.*def,'y.')
% plot(-1.2280232e-7*1000*strokedot.*span.*def.*defdot,'r.')
% plot(- 1.1480129e-6*1000.*span.*devddot,'g.')

%----
% width = 5.8;height = 2;
% ylabel('Force (N)')
% axis([cross(1) length(force) ylimit(1) ylimit(2)])
% set(gcf,'units','inches','papersize',[width height],'position',[1 9 width height],'paperpositionmode','auto')
% set(gca,'box','on','layer','top','Xtick',[])
% saveFolder = '/Users/charlesrichter/Desktop/CR Shared Drive/';
% saveas(gcf,[saveFolder,'test1'],'pdf')
%----


%% Angle Conventions and Angle of Attack
% To convert from my stroke convention to Sane/Dickinson, multiply by -1.
% phi_sd = -stroke; phidot_sd = -strokedot; phiddot_sd = -strokeddot;
% psi_sd = -def; psidot_sd = -defdot; psiddot_sd = -defddot;
% theta = dev; thetadot = devdot; thetaddot = devddot;

% defzjw = -1*def + pi/2;defzjwdot = -defdot;defzjwddot = -defddot;
% theta = def+pi/2; thetadot = defdot; thetaddot = defddot; % ZJW theta angle
% dev = zeros(size(dev));devdot = zeros(size(devdot));devddot = zeros(size(devddot));
% defrw = -1*(def-pi/2);defrwdot = -defdot;defrwddot = -defddot; % RJW is now the
% default angle convention.  This code was used to convert from ZJW to RJW

%% Angle of Attack (alpha)
% betaold = atan2(strokedot,devdot);
% alphaold = betaold-def;
% 
% beta_sd = atan2(phidot_sd,devdot);
% alpha_sd = psi_sd-beta_sd;

%% Translational Lift
% My Angles
% Clold = Clmax*sin(2*alphaold);
% Lold = .5*rho*chord*(1/3)*(outerspan^3-innerspan^3).*Clold.*(strokedot.^2+devdot.^2).*sin(betaold);
% % Sane/Dickinson Angles
% Cl_sd = 1.7*sin(2*alpha_sd);
% L_sd =
% .5*rho*chord*(1/3)*(outerspan^3-innerspan^3).*Cl.*sqrt(phidot_sd.^2+thetadot.^2).*sin(-beta_sd);

%% EQ Bad New - Trained on Exp 95
% THIS IS THE EQ BAD EQUATION USED IN THE THESIS!!
% eqbad = (0.33191672.*def.*sin(def)./(5.3005781 + def.*devddot.*cos(def).*sin(def) + 503.32355.*dev - 0.69517928.*def - stroke) - 0.0026070664.*devdot - 0.0042015365)./cos(0.16622099 + devdot);

% 16
% eqbad = -0.00174423*1 + -0.00377575*def + -7.86766e-5*devddot + 0.0368637*def.*def + -3.27239e-7*devddot.*devddot;
% eqbad = (1.91452*sin(def).*sin(2.52155*def) - 0.344079*def)./(96.7113 + devddot) - 0.00233856;
% eqbad = (1.83394*def.*sin(2.54819*def) - 0.335569*def)./(96.6771 + devddot) - 0.00224014;
% eqbad = (10.8758 + 13.7638*def.*def - 1.39296*def)./(397.721 + devddot + 16.9369*devdot + devdot.*devddot + 16.9369*def.*devdot) - 0.0289569;

eqbad = (14.9488 + 13.5087*def.*def - 1.38486.*def)./(407.258 + devddot + 17.8828.*devdot + cos(2.69241 - 42.9085./(13.5087.*def.*def) - 163.587.*def)) - 0.0381286;
% eqbad2 = (13.2404 + 12.9665*def.*def - 1.2865.*def)./(397.817 + devddot + 18.2332.*devdot + cos(54.7954*def.*def - def)) - 0.0346866;
% eqbad3 = (12.6842 + 13.6818*def.*def - 1.37269.*def)./(395.801 + devddot + sin(devddot) + 17.2074*devdot) - 0.0336;

% eqbad = (10.8758 + 13.7638*def.*def - 1.39296*def)./(397.721 + devddot + 16.9369*devdot + devdot.*devddot + 16.9369*def.*devdot) - 0.0289569;
% eqbad2 = (12.4958 + 13.753*def.*def - 1.40138*def)./(395.335 + devddot + sin(devddot) + 16.6449*devdot) - 0.033191;
% eqbad3 = (10.9152 + 14.4129*def.*def - 1.44433*def)./(396.623 + devdot + devddot) - 0.029252;


% % EQ's trained on 43
% eqbad = (15.17094*stroke + 0.83807027*strokedot.*strokedot - ...
% 1.3671925*strokedot.*cos(0.91728598*strokedot) - ...
% 0.44860578*strokedot)./(4753.5811 + 12.157823*strokedot.*strokedot) - ...
% 0.0066235247;

% 42
% eqbad =  (stroke.*stroke + stroke.*stroke.*dev + cos(devdot - 0.1436) - 0.0397754.*strokedot.*devdot).*(6.72859e-6.*strokedot + 0.000925379.*strokedot.*def - 0.000147078.*devddot - stroke.*stroke.*dev - 0.00312618.*devdot);
% eqbad2 = (0.000150184.*devddot + stroke.*stroke + cos(devdot - 0.202283) - 0.0456788.*strokedot.*devdot).*(1.09027e-5.*strokedot + 0.000914577.*strokedot.*def - 0.000150184.*devddot - stroke.*stroke.*dev - 0.00327404.*devdot);
% eqbad3 = (stroke.*stroke + cos(0.0885422 - devdot) - 0.511168.*def.*devdot).*(0.000900826.*strokedot.*def + 8.11487e-7.*strokedot.*strokedot - 0.000153206.*devddot - stroke.*stroke.*dev - 0.00314395.*devdot);
% eqbad4 = (0.000931815.*strokedot.*def - 0.000155455.*devddot - stroke.*stroke.*dev - 0.0028632.*devdot)./cos(-1.19131.*stroke - 0.0734576);
% eqbad5 = (0.000939599.*strokedot.*def - 0.000155704.*devddot - stroke.*stroke.*dev - 0.00296407.*devdot)./(cos(stroke) - 0.00240125.*devddot - 0.00240125.*strokedot);
% eqbad6 = (0.000964254.*strokedot.*def + 0.0427816.*stroke.*dev - 0.000155975.*devddot - stroke.*stroke.*dev - 0.00254397.*devdot)./(cos(stroke) - 0.00278222.*devddot - 0.00278222.*strokedot);
% eqbad7 = (0.000938313.*strokedot.*def - 0.000157845.*devddot - stroke.*stroke.*dev - 0.00316002.*devdot)./cos(stroke);
% eqbad8 = (0.000947451.*strokedot.*def + 0.0258625.*stroke.*dev - 0.000153212.*devddot - stroke.*stroke.*dev - 0.0026509.*devdot)./(cos(stroke - 0.0258625) - 0.0026509.*devddot - 0.0026509.*strokedot);
% eqbad = (0.000940495.*strokedot.*def - 0.000154229.*devddot - stroke.*stroke.*dev - 0.00298915.*devdot)./cos(1.08005.*stroke);

% 84
% eqbad =  (0.000553591*stroke - 0.011369*def.*sin(0.0883815 - 2.08476*def - 0.0883815*stroke.*def - def.*def - 0.431856*stroke) - 0.0019351)./(0.987596 + 0.223816*def);
% eqbad = (0.000606752*stroke - 0.0115711*def.*sin(0.0779974 - 2.18945*def - 0.793571*def.*def - 0.454391*stroke) - 0.00204384)./(1.04091 + 0.13929*def);
% eqbad = 0.000613254.*stroke - 0.000613178.*def - 0.0111796.*def.*sin(-2.21581.*def - 0.662951.*def.*def - 0.451466.*stroke) - 0.00199367;
% eqbad = (0.472924 + 0.252383.*stroke.*stroke - 0.0282963.*stroke).*(0.0282963.*stroke + 0.787534.*sin(def).*sin(2.50069.*def + 0.680791.*stroke - 0.0831549) + 0.0138744.*sin(2.50069.*def + 0.680791.*stroke - 0.0831549) - 0.252383.*stroke.*stroke)./(51.7283 + devddot);

% % EQ's trained on Experiment 94
% eqbad = 0.38975877./(27.224819 + defdot + devddot + 1.9677156./def + 6.8947473./(def.*def)) - 0.0029908489;

% eqbad = 0.0103059 + (def.*def - 0.0512405*def - 0.527118)./(46.8844 - devddot);
% % eqbad = (1.4266*def.*def - 0.0747445*def - 0.812839)./(63.1636 - devddot) + 0.016513*cos(devdot) - 0.00476817;
% % eqbad = 0.0103399 + (def.*def - 0.0604302*devdot - 0.0604302*def - 0.516703)./(45.4537 - devddot);
% % % eqbad = 0.0462103 + 0.0111978*cos(devdot) - 0.0597439*cos(def - 0.024308);
% % eqbad = 0.0386648 - 0.0386648./(1.92731 + cos(3.10529 + 1.33415*def));
% 
% eqbad = (1.26907*def.*def.*cos(devdot) - 0.0645029*def - 0.676606)./(56.7135 - devddot) + 0.0143676*cos(devdot) - 0.00351427;
% eqbad = 1.10967*cos(devdot).*sin(0.0224027 + (def.*def - 0.0211669*devdot - 0.0530431*def - 0.545035)./(49.3636 - 0.861782*devddot)) - 0.0137155;
% eqbad = (0.35772*def.*def.*devdot + def.*def - 0.571619*cos(devdot) - def.*def.*def.*devdot - 0.0456542*def)./(48.9517 - devddot) + 0.0297071*cos(devdot) - 0.0187454;
% eqbad = (1.1581*def.*def.*cos(devdot) + 1.1581*def.*devdot.*devdot - 0.124974*def.*devdot - 0.675676*cos(devdot) - 0.0639685*def)./(54.308 - devddot) + 0.0378883*cos(devdot) - 0.0263563;
% 
% eqbad = 1.10967*cos(devdot).*sin(0.0224027 + (def.*def - 0.0211669*devdot - 0.0530431*def - 0.545035)./(49.3636 - 0.861782*devddot)) - 0.0137155;
% eqbad = 0.0120687 + (0.205368*def.*def - 0.0120687*def - 0.121603)./(9.20141 - 0.149044*devddot - devdot);
% eqbad = 0.00494643 - 0.00790684*cos(0.0828926 - 3.10557*def);
% 
% % Report 4
% eqbad = (1.06712.*def.*def.*cos(devdot) - 0.0731058.*def.*devdot - 0.54871.*cos(devdot) - 0.0528022.*def)./(49.1837 - 1.05787.*devddot) + 0.0230924.*cos(devdot) - 0.0128708;
% eqbad = 1.10967.*cos(devdot).*sin(0.0224027 + (def.*def - 0.0211669.*devdot - 0.0530431.*def - 0.545035)./(49.3636 - 0.861782.*devddot)) - 0.0137155;
% eqbad = (1.02156.*def.*def.*cos(devdot) - 0.049263.*def.*cos(devdot) - 0.0944365.*def.*devdot.*cos(devdot) - 0.549786.*cos(devdot))./(46.8129 - 0.917294.*devddot) + 0.023672.*cos(devdot) - 0.0129301;
% eqbad = 1.1219.*cos(devdot).*sin(0.0223679 + (0.974298.*def.*def - 0.0199563.*devdot - 0.0491337.*def - 0.543857)./(49.8849 - 0.913784.*devddot)) - 0.0138736;
% eqbad = (1.26907.*def.*def.*cos(devdot) - 0.0645029.*def - 0.676606)./(56.7135 - devddot) + 0.0143676.*cos(devdot) - 0.00351427;
% eqbad = (1.29713.*def.*def - 0.0661413.*def - 0.689287)./(57.6032 - devddot) + 0.0177224.*cos(devdot) - 0.0068506;
% 
% % Report 3
% eqbad = 1.10967.*cos(devdot).*sin(0.0224027 + (def.*def - 0.0211669.*devdot - 0.0530431.*def - 0.545035)./(49.3636 - 0.861782.*devddot)) - 0.0137155;
% eqbad = (1.26907.*def.*def.*cos(devdot) - 0.0645029.*def - 0.676606)./(56.7135 - devddot) + 0.0143676.*cos(devdot) - 0.00351427;
% eqbad =  (1.29713.*def.*def - 0.0661413.*def - 0.689287)./(57.6032 - devddot) + 0.0177224.*cos(devdot) - 0.0068506;
% eqbad = 0.0120687 + (0.205368.*def.*def - 0.0120687.*def - 0.121603)./(9.20141 - 0.149044.*devddot - devdot);
% eqbad = 0.0103059 + (def.*def - 0.0512405.*def - 0.527118)./(46.8844 - devddot);
% 
% % Report 2
% eqbad = (1.23103*def.*def.*cos(devdot) - 0.0626916*def - 0.644645)./(54.9504 - devddot) + 0.0146582*cos(devdot) - 0.00399526;
% % eqbad2 = (1.22827*def.*def - 0.0620891*def - 0.658202)./(54.9212 - devddot) + 0.0149733*cos(devdot) - 0.00408114;
% % eqbad2 = 0.0120687 + (0.205368*def.*def - 0.0120687*def - 0.121603)./(9.20141 - 0.149044*devddot - devdot);
% 
% % eqbad2 = 0.0120687 + (0.205368*def.*def - 0.0120687*def - 0.121603)./(9.20141 - 0.149044*devddot - devdot);
% % eqbad3 =  0.0103399 + (def.*def - 0.0604302*devdot - 0.0604302*def - 0.516703)./(45.4537 - devddot);
% 
% % eqbad2 =  0.0321381*sin(def).*sin(def) + 0.0351362*def.*def.*devdot - 0.00554033*def.*devdot - 0.00492998*devdot - 0.00140979*def - 0.00246132;
% % eqbad3 = 0.030665.*def.*def.*devdot + 0.030665.*def.*sin(def) - 0.00508014.*def.*devdot - 0.00432505.*devdot - 0.00145131.*def - 0.00242072;
% 
% % eqbad =  1.1219.*cos(devdot).*sin(0.0223679 + (0.974298.*def.*def - 0.0199563.*devdot - 0.0491337.*def - 0.543857)./(49.8849 - 0.913784.*devddot)) - 0.0138736;
% % eqbad2 = (1.26907.*def.*def.*cos(devdot) - 0.0645029.*def - 0.676606)./(56.7135 - devddot) + 0.0143676.*cos(devdot) - 0.00351427;
% % eqbad3 = 0.0123666.*cos(devdot) + (1.32711.*def.*def - 0.0664403.*def - 0.808785)./(60.1618 - devddot);

% % Report 1
% eqbad = 0.0103016 + (0.834639*def.*def - 0.0448558*def - 0.435238)./(38.3609 - 0.760572*devddot - devdot);
% eqbad = 0.0321381*sin(def).*sin(def) + 0.0351362*def.*def.*devdot - 0.00554033*def.*devdot - 0.00492998*devdot - 0.00140979*def - 0.00246132;
% eqbad = 0.030665*def.*def.*devdot + 0.030665*def.*sin(def) - 0.00508014*def.*devdot - 0.00432505*devdot - 0.00145131*def - 0.00242072;
% eqbad = 0.0114767 + (0.851299*def.*def - 0.0465245*def - 0.482987)./(38.2765 - 0.658899*devddot);

% figure(1),clf,hold on
% set(gcf,'defaultlinelinewidth',1)
% plot(force,'k')
% plot(eqbad,'b')
% ylimit = axbuf.*[min([force;(eqbad)]) max([force;(eqbad)])];
% ylim(ylimit)
% xlim([1 length(force)])
% set(gca,'box','on','layer','top','Xtick',[])
% set(gcf,'units','inches','papersize',[2 1],'position',[1 9 2 1],'paperpositionmode','auto')
% saveFolder = '/Users/charlesrichter/Desktop/CR Shared Drive/Thesis/figures/';
% saveas(gcf,[saveFolder,'eqnew_bad_',num2str(toplot(bla))],'pdf')

%% EQ Bad
% eqbad = 0.00044409049 + 1.2689823e-5*stroke.*strokeddot + 0.29627392*dev - 1.9190562e-5*def.*strokeddot - 1.6229924e-6*strokeddot - 0.0027215057*def;

%% EQ Bad 2
% eqbad2 = 0.022154899 + (0.0047265203./cos(devdot) - 0.025963018*cos(0.10153064 + def))./cos(stroke.*stroke - 0.10847107)...
%     - 0.0038038024*stroke.*stroke.*devdot - 0.0038038024*stroke.*stroke.*def;

% beta0 = [0.022154899 0.0047265203 0.025963018 0.10153064 0.10847107 0.0038038024 0.0038038024];
% beta = nlinfit(eqdata,force,@fn_eqbad,beta0);
% forcepredict = fn_eqbad(beta,eqdata);

%% EQ 1
% Leq1 = 1.4488915e-8*1000*span.*Ycm.*defdot.*defdot.*cos(def.*def - def).*cos(def.*def - def); % -- ORIGINAL COEFFS
% Leq1 = .8e-8*1000*span.*Ycm.*defdot.*defdot.*cos(def.*def - def).*cos(def.*def - def);
% Leq2 = 9.7118491e-9*strokedot*1000.*span.*Xcm;
% Leq3 = 1.2148172e-7*strokedot*1000.*span.*Xcm.*def; % -- ORIGINAL COEFFS
% Leq3 = 1e-7*strokedot*1000.*span.*Xcm.*def;
% Leq4 = -1.2165592e-6*span*1000.*devddot;

% figure(1),clf,hold on
% ylimit = axbuf*[min(force) max(force)];
% cross = find([0;diff(sign(strokedot))]);
% fill([0 cross(1) cross(1) 0],[-1 -1 1 1],shadecolor,'edgecolor',shadecolor)
% fill([cross(2) cross(3) cross(3) cross(2)],[-1 -1 1 1],shadecolor,'edgecolor',shadecolor)
% set(gcf,'defaultlinelinewidth',1)
% plot(domain,force,'k')
% plot(domain,inertialdev,'c')
% plot(domain,inertial,'b')
% plot(domain,inertial+inertialdev,'r')

% width = 5.8;height = 2;
% ylabel('Force (N)')
% axis([1 length(force) ylimit(1) ylimit(2)])
% set(gcf,'units','inches','papersize',[width height],'position',[-8 9 width height],'paperpositionmode','auto')
% set(gca,'box','on','layer','top','Xtick',[])
% saveFolder = '\\zeus\shared\Charlie Richter\Thesis\figures\';
% saveas(gcf,[saveFolder,'inertial1'],'pdf')

%% EQ new to use: 2 (1),4 (2),5 (3),7 (4),8 (5),9 (6)
%% EQ New 1
% 1.6719996e-8*span*Ycm*defdot*defdot*cos(def)*cos(def) + 1.1938569e-7*strokedot*span*Xcm*def - 1.1938569e-7*strokedot*span*def*devddot - 1.1938569e-7*strokedot*span*def*defdot - 6.8213382e-7*span*devddot
% figure(1),clf,hold on,set(gcf,'defaultlinelinewidth',2),plot(force,'k'),plot(1.6719996e-8.*1000.*span.*Ycm.*defdot.*defdot.*cos(def).*cos(def),'r')
%% EQ New 2 = 1
% 1.5723863e-8*span*Ycm*defdot*defdot*cos(def) + 1.2146853e-7*strokedot*span*Xcm*def - 1.2146853e-7*strokedot*span*def*devddot - 1.2146853e-7*strokedot*span*def*defdot - 7.0700918e-7*span*devddot
% figure(1),hold on,set(gcf,'defaultlinelinewidth',2),plot(force,'k'),plot(1.5723863e-8.*1000.*span.*Ycm.*defdot.*defdot.*cos(def),'b')
%% EQ New 3
% 1.5422589e-8*span*Ycm*defdot*defdot*cos(def)*cos(def) + 1.223172e-7*strokedot*span*Xcm*def - 1.223172e-7*strokedot*span*def*defdot - 1.0904402e-6*span*devddot
%% EQ New 4 = 2
% 1.419276e-8*span*Ycm*defdot*defdot*cos(def) +
% 1.2280232e-7*strokedot*span*Xcm*def - 1.2280232e-7*strokedot*span*def*defdot - 1.1480129e-6*span*devddot
% eq1 = 1.419276e-8.*1000..*span.*Ycm.*defdot.*defdot.*cos(def);
% eq2 = 1.2280232e-7.*strokedot.*1000.*span.*Xcm.*def;
% eq3 = -1.2280232e-7.*strokedot.*1000.*span.*def.*defdot;
% eq4 = -1.1480129e-6.*1000.*span.*devddot;
% 
% figure(1),clf,hold on
% plot(force)
% plot(eq1+eq2+eq3+eq4)
% plot(eq2)
% plot(Lrjw)

%% EQ New 5 = 3
% 1.4180717e-8*span*Ycm*defdot*defdot*cos(def - def*def) + 1.2409225e-7*strokedot*span*Xcm*def - 1.0660193e-6*span*devddot
% figure(1),clf,hold on,set(gcf,'defaultlinelinewidth',2),plot(force,'k'),plot(- 1.0660193e-6*1000*span.*devddot,'b')
%% EQ New 6
% 1.4679481e-8*span*Ycm*defdot*defdot*cos(def - 0.51079893) + 1.2142625e-7*strokedot*span*Xcm*def - 1.1698893e-6*span*devddot
% figure(1),clf,hold on,set(gcf,'defaultlinelinewidth',2),plot(force,'k'),plot(1.4679481e-8*1000*span.*Ycm.*defdot.*defdot.*cos(def - 0.51079893),'b')
%% EQ New 7 = 4
% 1.4526972e-8*span*Ycm*defdot*defdot*cos(def) + 1.2489714e-7*strokedot*span*Xcm*def - 1.0872654e-6*span*devddot
% eq1 = 1.4526972e-8.*1000.*span.*Ycm.*defdot.*defdot.*cos(def);
% eq2 = 1.2489714e-7.*1000.*strokedot.*span.*Xcm.*def;
% eq3 = -1.0872654e-6.*1000.*span.*devddot;
% 
% figure(1),clf,hold on
% plot(force)
% plot(eq1+eq2+eq3)
%% EQ New 8 = 5
% 1.2524723e-7*strokedot*span*Xcm*def - 5.8678208e-11*span*lowerchord*lowerchord*defdot*defdot - 1.1539232e-6*span*devddot
%% EQ New 9 = 6
% 1.3552796e-8*span*Ycm*defdot*defdot + 1.2505514e-7*strokedot*span*Xcm*def - 1.1357268e-6*span*devddot
%% EQ New 10
% 1.320138e-8*span*Ycm*defdot*defdot + 1.2529193e-7*strokedot*span*Xcm*def - 9.836931e-5*devddot
%% EQ New 11
% 1.2343477e-8*Iyx*defdot*defdot + 9.7913826e-6*strokedot*span*def - 1.1256892e-6*span*devddot
% figure(1),clf,hold on,set(gcf,'defaultlinelinewidth',2),plot(force,'k'),plot(9.7913826e-6*1000*strokedot.*span.*def,'b')
%% EQ New 12
% 1.2400019e-8*Ixy*defdot*defdot + 9.7674611e-6*strokedot*span*def - 0.00010300481*devddot
%% EQ New 13
% 9.6821368e-6*strokedot*span*def - 1.02371e-5*defdot*defdot - 0.00010280222*devddot
%% EQ New 14
% 1.0339508e-5*strokedot*span*def + 1.4668144e-6*Ycm*defdot*defdot
%% EQ New 15
% 1.0226605e-5*strokedot*span*def - 1.2614417e-5*defdot*defdot
%% EQ New 16
% 1.0738358e-7*strokedot*def*Iyy - 0.00015105731*devddot

%% RJW Calculation
Clmax = 1.9;
beta_trans = 1;
beta_rot = 1;

% % Expressions for wy,wz taken from RJW equation 2.8, page 200
wx = defdot-strokedot.*sin(dev);
wxdot = defddot - (strokeddot.*sin(dev)+strokedot.*devdot.*cos(dev));
wz = devdot.*cos(def) + strokedot.*sin(def).*cos(dev); % Vo = r_i*omega_z
wy = -strokedot.*cos(def).*cos(dev) + devdot.*sin(def); % Wo = r_i*omega_y
wydot = -strokeddot.*cos(def).*cos(dev)+strokedot.*defdot.*sin(def).*cos(dev)+...
    strokedot.*devdot.*cos(def).*sin(dev)+devddot.*sin(def)+devdot.*defdot.*cos(def);
wh = sqrt(wz.^2+wy.^2); % angular velocity of the wing hinge

alpha = atan2(-wy,wz); % angle of attack relative to instantaneous velocity
beta = alpha+def; % angle between velocity vector and vertical reference
Cl = Clmax*sin(2*alpha); % coefficient of lift
Lrjw = .5*rho.*chord*(1/3).*(outerspan.^3-innerspan.^3).*Cl.*wh.^2.*sin(beta); % kg/m^3 * rad^2/s^2 * m * m^3 = kg*m/s^2

Z01 = beta_trans*-lambdaz*.5.*(wydot-wx.*wz).*-(outerspan.^2-innerspan.^2);
Z02 = beta_rot*-lambdazw.*wxdot.*(outerspan-innerspan);
Z0spanwise = Z01+Z02;
AMrjw = Z0spanwise.*-sin(def);
RJWtot = Lrjw+AMrjw+inertial+inertialdev;

% figure(2),clf,hold on
% ylimit = axbuf*[min([force;RJWtot]) max([force;RJWtot])];
% cross = find([0;diff(sign(strokedot))]);
% fill([0 cross(1) cross(1) 0],[-1 -1 1 1],shadecolor,'edgecolor',shadecolor)
% fill([cross(2) cross(3) cross(3) cross(2)],[-1 -1 1 1],shadecolor,'edgecolor',shadecolor)
% 
% set(gcf,'defaultlinelinewidth',1)
% plot(domain,force,'k')
% plot(domain,Lrjw+AMrjw+inertial+inertialdev,'r')
% plot(domain,Lrjw,'b')
% plot(domain,Z01.*-sin(def),'g')
% plot(domain,Z02.*-sin(def),'m')
% 
% width = 5.8;height = 2;
% ylabel('Force (N)')
% axis([1 length(force) ylimit(1) ylimit(2)])
% set(gcf,'units','inches','papersize',[width height],'position',[-8 6 width height],'paperpositionmode','auto')
% set(gca,'box','on','layer','top','Xtick',[])
% saveFolder = '\\zeus\shared\Charlie Richter\Thesis\figures\';
% saveas(gcf,[saveFolder,'forces_rw_am'],'pdf')

%% Sane/Dickinson Rotational Lift Using RJW Coordinates
% Frot2 = Crot*rho*chord^2*.5*(outerspan^2-innerspan^2).*wh.*-wx.*sin(def);

%% ZJW Calculation
% gamma = -2*Ct*a*RRR*wz.*-wy./sqrt(wz.^2+wy.^2)+2*Cr*a^2*defdot;
% Fv = rho*a*(A-B*(wz.^2-wy.^2)./(wz.^2+wy.^2)).*sqrt(wz.^2+wy.^2)*RRR;
% Fxprime = m22*defdot*RRR*-wy - rho*RRR*-wy*gamma - Fv.*RRR.*wz;
% RRR represents multiplication by r -> stick a dr on the end and integrate
% Expand to:
% Fxprime = m22*defdot*RRR*-wy -...
%     rho*-wy.*(-2*Ct*a*RRR*RRR*wz.*-wy./sqrt(wz.^2+wy.^2)+2*Cr*a^2*RRR*defdot) -...
%     rho*a*(A-B*(wz.^2-wy.^2)./(wz.^2+wy.^2)).*sqrt(wz.^2+wy.^2)*RRR*RRR.*
%     wz;

Ct = 1.69;
Cr = 1;
m22 = .75*lambdaz;

R2 = (outerspan.^2-innerspan.^2)/2;
R3 = (outerspan.^3-innerspan.^3)/3;
Fxprime = m22.*defdot.*R2.*-wy -...
    rho*-wy.*(-2*Ct.*a.*R3.*wz.*-wy./sqrt(wz.^2+wy.^2)+2*Cr.*a.^2.*R2.*defdot);
Fyprime = -m22.*R2.*(-wydot+wx.*wz) +...
    rho*wz.*(-2*Ct.*a.*R3.*wz.*-wy./sqrt(wz.^2+wy.^2)+2*Cr.*a.^2.*R2.*defdot);
Lzjw = Fxprime.*cos(def)-Fyprime.*sin(def);
ZJWtot = Lzjw+inertial+inertialdev;

% Examine individual components of ZJW model
% Added mass term 1 - NUMERICALLY IDENTICAL TO RJW TERM: Z01.*-SIN(DEF)
Fxprime1 = zeros(length(Fxprime),1);
Fyprime1 = -m22.*R2.*(-wydot+wx.*wz);
L1 = Fxprime1.*cos(def)-Fyprime1.*sin(def);

% Added mass term 2
Fxprime2 = m22.*defdot.*-wy.*R2;
Fyprime2 = zeros(length(Fxprime),1);
L2 = Fxprime2.*cos(def)-Fyprime2.*sin(def);

% Circulation term 1 - NUMERICALLY IDENTICAL TO RJW LIFT TERM: Lrjw 
Fxprime3 = -rho*-wy.*R3*-2*Ct.*a.*wz.*-wy./sqrt(wz.^2+wy.^2); % First half of Circulation term
Fyprime3 = rho*wz.*R3*-2*Ct.*a.*wz.*-wy./sqrt(wz.^2+wy.^2); % First half of Circulation term
L3 = Fxprime3.*cos(def)-Fyprime3.*sin(def);

% Circulation term 2 - SIMILAR TO SANE/DICKINSON ROTATIONAL LIFT
Fxprime4 = -rho*-wy.*R2*2*Cr.*a.^2.*defdot; % Second half of Circulation term
Fyprime4 = rho*wz.*R2*2*Cr.*a.^2.*defdot; % Second half of Circulation term
L4 = Fxprime4.*cos(def)-Fyprime4.*sin(def);

% % Viscous term - Almost zero contribution
% Fxprime5 = -rho*a.*(A-B*(wz.^2-wy.^2)./(wz.^2+wy.^2)).*sqrt(wz.^2+wy.^2).*wz.*R3;
% Fyprime5 = -rho*a.*(A-B*(wz.^2-wy.^2)./(wz.^2+wy.^2)).*sqrt(wz.^2+wy.^2).*-wy.*R3;
% L5 = Fxprime5.*cos(def)-Fyprime5.*sin(def);

% figure(3),clf,hold on
% ylimit = axbuf*[min([force;ZJWtot]) max([force;ZJWtot])];
% cross = find([0;diff(sign(strokedot))]);
% fill([0 cross(1) cross(1) 0],[-1 -1 1 1],shadecolor,'edgecolor',shadecolor)
% fill([cross(2) cross(3) cross(3) cross(2)],[-1 -1 1 1],shadecolor,'edgecolor',shadecolor)

% set(gcf,'defaultlinelinewidth',1)
% plot(domain,force,'k')
% plot(domain,ZJWtot,'r')
% plot(Lzjw+inertial+inertialdev,'r')
% plot(L3,'b')
% plot(L4,'c')
% plot(L1,'g')
% plot(L2,'m')
% 
% width = 5.8;height = 2;
% ylabel('Force (N)')
% axis([1 length(force) ylimit(1) ylimit(2)])
% set(gcf,'units','inches','papersize',[width height],'position',[-8 3 width height],'paperpositionmode','auto')
% set(gca,'box','on','layer','top','Xtick',[])
% saveFolder = '\\zeus\shared\Charlie Richter\Thesis\figures\';
% saveas(gcf,[saveFolder,'forces_jw_am'],'pdf')

%% EQ New 7 = 4
% 1.4526972e-8*span*Ycm*defdot*defdot*cos(def) + 1.2489714e-7*strokedot*span*Xcm*def - 1.0872654e-6*span*devddot

eq1 = 1.4526972e-8.*1000.*span.*Ycm.*defdot.*defdot.*cos(def);
C2 = 1.2489714e-7;
eq2 = C2.*1000.*strokedot.*span.*Xcm.*def;
eq3 = -1.0872654e-6.*1000.*span.*devddot;

% inertial = -yh.*chord.*thickness.*span.*rho720.*(defddot.*sin(def) + defdot.^2.*cos(def));
% inertialdev = -xh.*span.*chord.*thickness.*rho720.*(devddot.*cos(dev) - devdot.^2.*sin(dev));

% figure(2),clf,hold on
% set(gcf,'defaultlinelinewidth',1)
% plot(force,'k')
% plot(eq1+eq2+eq3,'r')
% ylimit = axbuf*[min([force;(eq1+eq2+eq3)]) max([force;(eq1+eq2+eq3)])];
% ylim(ylimit)
% xlim([1 length(force)])
% set(gca,'box','on','layer','top','Xtick',[])
% set(gcf,'units','inches','papersize',[2 1],'position',[3 9 2 1],'paperpositionmode','auto')
% saveFolder = '/Users/charlesrichter/Desktop/CR Shared Drive/Thesis/figures/';
% saveas(gcf,[saveFolder,'eqnew_good_',num2str(toplot(bla))],'pdf')


%----
if bla == 5
    
fontsz = 6;
widthcm = 8.7;
heightcm = 6;
roundmult = .001;
plotfactor = 1000; % Multiply force by 1000 then explain that all forces are "e-3"
plottitles = {'Training Data','Validation Experiment 1','Validation Experiment 2','Validation Experiment 3'};

% Buffers in CM
lbuf = 1;
rbuf = 1;
hbuf = .6;
tbuf = .2;
bbuf = 1;
vbuf = .6;

nhorz = 1;
nvert = 2;

axwidth = (widthcm - lbuf - rbuf - (nhorz-1)*hbuf)/nhorz;
horzgrid = linspace(lbuf,lbuf + (nhorz-1)*(axwidth + hbuf),nhorz);
axheight = (heightcm - tbuf - bbuf - (nvert-1)*vbuf)/nvert;
vertgrid = fliplr(linspace(bbuf,bbuf + (nvert-1)*(axheight+vbuf),nvert));

counter = 1;
for nn = 1:length(vertgrid)
    for mm = 1:length(horzgrid)
        plotgrid1{counter} = [horzgrid(mm) vertgrid(nn)]; %#ok<SAGROW>
        counter = counter+1;
    end
end
           
% fontsz = 6;
% widthcm = 8.7;
% heightcm = 5.6;
% plotfactor = 1000; % Multiply force by 1000 then explain that all forces are "e-3"
% 
% firstbuf = 0.1;%*17.8/8.7;
% horzbuf = .04;%*17.8/8.7;
% vertbuf = .1;
% topvertbuf = .06;
% bottomvertbuf = .04;
% axwidth = 1-firstbuf-horzbuf;
% axheight = (1-bottomvertbuf-topvertbuf-vertbuf)/2;
% 
figure(1)
set(gcf,'units','centimeters','paperunits','centimeters','papersize',[widthcm heightcm],...
    'paperposition',[0 0 widthcm heightcm])
set(gcf,'defaultlinelinewidth',.5,'visible','off')

% subplot(2,1,1),hold on
axes('units','centimeters','position',[plotgrid1{1} axwidth axheight]),hold on
cross1 = find([0;diff(sign(strokedot))]);

normstroke = [1;cross1;length(force)];
xticklab = {'0','.5','1.0','1.5','2'};

fill([0 cross1(1) cross1(1) 0],10000*[-1 -1 1 1],shadecolor,'edgecolor',shadecolor)
fill([cross1(2) cross1(3) cross1(3) cross1(2)],10000*[-1 -1 1 1],shadecolor,'edgecolor',shadecolor)
plot(force*plotfactor,'k')
plot(eq2*plotfactor,'b')
plot(Lrjw*plotfactor,'r')
ylimit = axbuf*[min([force*plotfactor;(eq1+eq3)*plotfactor]) max([force*plotfactor;(eq1+eq3)*plotfactor])];
ylim(ylimit)
xlim([1 length(force)])
% set(gca,'units','normalized','position',[firstbuf bottomvertbuf+vertbuf+axheight axwidth axheight])
% axis([1 length(force) ylimit(1) ylimit(2)])
set(gca,'box','on','layer','top','Xtick',normstroke,'XTickLabel',xticklab,'FontUnits','points','FontSize',fontsz)
ylabel('Force (10^-^3 N)','FontSize',fontsz)
text(axwidth-.24,axheight-.1,'A','units','centimeters',...
    'HorizontalAlignment','right','VerticalAlignment','top',...
    'fontsize',10)

% subplot(2,1,2),hold on
axes('units','centimeters','position',[plotgrid1{2} axwidth axheight]),hold on
cross1 = find([0;diff(sign(strokedot))]);
fill([0 cross1(1) cross1(1) 0],10000*[-1 -1 1 1],shadecolor,'edgecolor',shadecolor)
fill([cross1(2) cross1(3) cross1(3) cross1(2)],10000*[-1 -1 1 1],shadecolor,'edgecolor',shadecolor)
plot(force*plotfactor,'k')
plot((eq1+eq3)*plotfactor,'b')
plot((inertial+inertialdev)*plotfactor,'r')
ylimit = axbuf*[min([force*plotfactor;(eq1+eq3)*plotfactor]) max([force*plotfactor;(eq1+eq3)*plotfactor])];
ylim(ylimit)
xlim([1 length(force)])
ylabel('Force (10^-^3 N)','FontSize',fontsz)
% axis([1 length(force) ylimit(1) ylimit(2)])
% set(gca,'units','normalized','position',[firstbuf bottomvertbuf axwidth axheight])
set(gca,'box','on','layer','top','Xtick',normstroke,'XTickLabel',xticklab,'FontUnits','points','FontSize',fontsz)
xlabel('Stroke Cycle')
text(axwidth-.24,axheight-.1,'B','units','centimeters',...
    'HorizontalAlignment','right','VerticalAlignment','top',...
    'fontsize',10)

% set(gcf,'units','inches','papersize',[6 6],'paperpositionmode','auto')
##saveFolder = '/Users/charlesrichter/Desktop/PNAS11_Richter/figures/';
saveFolder = '/home/crichter/ccsl/NSR20_Richter/paper/figures/';
saveas(gcf,[saveFolder,'eq_comparison_PNAS'],'epsc')
end

%% Plot for PNAS
if bla < 5
txbuf = .2;
fontsz = 6;
widthcm2 = 17.8;
heightcm2 = 6;
roundmult = .001;
plotfactor = 1000; % Multiply force by 1000 then explain that all forces are "e-3"
plottitles = {'Training Data','Validation Experiment 1','Validation Experiment 2','Validation Experiment 3'};
plotlabels = {'A','B','C','D'};
xticklab = {'0','.5','1.0','1.5','2'};

% firstbuf = .08;
% horzbuf = .04;
% vertbuf = .1;
% topvertbuf = .08;
% bottomvertbuf = .02;
naxes = 4;

% axwidth = (1-(firstbuf + naxes*horzbuf))/naxes;
% plotgrid = linspace(firstbuf,firstbuf+3*(axwidth+horzbuf),4);
% axheight = (1-bottomvertbuf-topvertbuf-vertbuf)/2;
% vertgrid = linspace(bottomvertbuf,bottomvertbuf+vertbuf+axheight,2);

% Buffers in CM
lbuf2 = 1;
rbuf2 = .2;
hbuf2 = .6;
tbuf2 = .2;
bbuf2 = 1;
vbuf2 = .6;

nhorz = 2;
nvert = 2;

axwidth2 = (widthcm2 - lbuf2 - rbuf2 - (nhorz-1)*hbuf2)/nhorz;
horzgrid2 = linspace(lbuf2,lbuf2 + (nhorz-1)*(axwidth2 + hbuf2),nhorz);
axheight2 = (heightcm2 - tbuf2 - bbuf2 - (nvert-1)*vbuf2)/nvert;
vertgrid2 = fliplr(linspace(bbuf2,bbuf2 + (nvert-1)*(axheight2+vbuf2),nvert));

counter = 1;
for nn = 1:length(vertgrid2)
    for mm = 1:length(horzgrid2)
        plotgrid{counter} = [horzgrid2(mm) vertgrid2(nn)]; %#ok<SAGROW>
        counter = counter+1;
    end
end
        
figure(3)
set(gcf,'units','centimeters','paperunits','centimeters','papersize',[widthcm2 heightcm2],...
    'paperposition',[0 0 widthcm2 heightcm2])
set(gcf,'defaultlinelinewidth',.5,'visible','off')

axes('units','centimeters','position',[plotgrid{bla} axwidth2 axheight2]),hold on

cross1 = find([0;diff(sign(strokedot))]);
if cross1(1)>50
    fill([0 cross1(1) cross1(1) 0],10000*[-1 -1 1 1],shadecolor,'edgecolor',shadecolor)
    fill([cross1(2) cross1(3) cross1(3) cross1(2)],10000*[-1 -1 1 1],shadecolor,'edgecolor',shadecolor)
else
    fill([0 cross1(2) cross1(2) 0],10000*[-1 -1 1 1],shadecolor,'edgecolor',shadecolor)
    fill([cross1(3) cross1(4) cross1(4) cross1(3)],10000*[-1 -1 1 1],shadecolor,'edgecolor',shadecolor)
end

normstroke = [1;cross1(length(cross1)-2:end);length(force)];

plot(force * plotfactor,'k')
plot(eqbad * plotfactor,'r')
plot((eq1+eq2+eq3) * plotfactor,'b')

text(axwidth2-.24,axheight2-.1,plotlabels{bla},'units','centimeters',...
    'HorizontalAlignment','right','VerticalAlignment','top',...
    'fontsize',10)

ylimit = axbuf*[min([force*plotfactor;(eqbad)*plotfactor;(eq1+eq2+eq3)*plotfactor])...
    max([force*plotfactor;(eqbad)*plotfactor;(eq1+eq2+eq3)*plotfactor])];
ylim(ylimit)
xlim([1 length(force)])
set(gca,'box','on','layer','top','Xtick',normstroke,'XTickLabel',xticklab,'FontUnits','points','FontSize',fontsz)
% title(plottitles(bla),'FontSize',fontsz)
if bla==1 || bla==3,ylabel('Force (10^-^3 N)'),end
if bla==3 || bla==4,xlabel('Stroke Cycle'),end
end


% % subplot(2,naxes,bla),hold on
% subplot(naxes,2,2*bla-1),hold on
% cross1 = find([0;diff(sign(strokedot))]);
% if cross1(1)>50
%     fill([0 cross1(1) cross1(1) 0],10000*[-1 -1 1 1],shadecolor,'edgecolor',shadecolor)
%     fill([cross1(2) cross1(3) cross1(3) cross1(2)],10000*[-1 -1 1 1],shadecolor,'edgecolor',shadecolor)
% else
%     fill([0 cross1(2) cross1(2) 0],10000*[-1 -1 1 1],shadecolor,'edgecolor',shadecolor)
%     fill([cross1(3) cross1(4) cross1(4) cross1(3)],10000*[-1 -1 1 1],shadecolor,'edgecolor',shadecolor)
% end
% plot(force * plotfactor,'k')
% plot(eqbad * plotfactor,'b')
% ylimit = axbuf*[min([force*plotfactor;(eqbad)*plotfactor]) max([force*plotfactor;(eqbad)*plotfactor])];
% % ylimit(1) = floor(ylimit(1)/roundmult)*roundmult;
% % ylimit(2) = ceil(ylimit(2)/roundmult)*roundmult;
% % ydiff = (ylimit(2) - ylimit(1))/roundmult;
% % Yticks = ylimit(1):roundmult:ylimit(2);
% % Yticks = linspace(ylimit(1),ylimit(2),8);
% ylim(ylimit)
% xlim([1 length(force)])
% % title(plottitles(bla),'FontSize',fontsz)
% 
% % set(gca,'position',[plotgrid(bla) vertgrid(2) axwidth axheight])
% set(gca,'units','centimeters','position',[plotgrid(1) vertgrid(naxes+1-bla) axwidth axheight]);
% set(gca,'box','on','layer','top','Xtick',[],'FontUnits','points','FontSize',fontsz)
% % set(gca,'ActivePositionProperty','position')
% % if bla==1,ylabel({'Single Training';'Force (10^-^3 N)'}),end
% 
% figure(3)
% % subplot(2,naxes,naxes+bla),hold on
% subplot(naxes,2,2*bla),hold on
% if cross1(1)>50
%     fill([0 cross1(1) cross1(1) 0],10000*[-1 -1 1 1],shadecolor,'edgecolor',shadecolor)
%     fill([cross1(2) cross1(3) cross1(3) cross1(2)],10000*[-1 -1 1 1],shadecolor,'edgecolor',shadecolor)
% else
%     fill([0 cross1(2) cross1(2) 0],10000*[-1 -1 1 1],shadecolor,'edgecolor',shadecolor)
%     fill([cross1(3) cross1(4) cross1(4) cross1(3)],10000*[-1 -1 1 1],shadecolor,'edgecolor',shadecolor)
% end
% plot(force * plotfactor,'k')
% plot((eq1+eq2+eq3) * plotfactor,'r')
% ylimit = axbuf*[min([force*plotfactor;(eq1+eq2+eq3)*plotfactor]) max([force*plotfactor;(eq1+eq2+eq3)*plotfactor])];
% % ylimit(1) = floor(ylimit(1)/roundmult)*roundmult;
% % ylimit(2) = ceil(ylimit(2)/roundmult)*roundmult;
% ylim(ylimit)
% xlim([1 length(force)])
% % set(gca,'position',[plotgrid(bla) vertgrid(1) axwidth axheight])
% set(gca,'units','centimeters','position',[plotgrid(2) vertgrid(naxes+1-bla) axwidth axheight],'units','centimeters')
% set(gca,'box','on','layer','top','Xtick',[],'FontUnits','points','FontSize',fontsz)
% % set(gca,'YTickLabel', num2str(get(gca,'YTick')','%g'))
% % if bla==1,ylabel({'Multiple Training';'Force (10^-^3 N)'}),end
% end

if bla == 4
##saveFolder = '/Users/charlesrichter/Desktop/PNAS11_Richter/figures/';
saveFolder = '/home/crichter/ccsl/NSR20_Richter/paper/figures/';
saveas(gcf,[saveFolder,'eq_goodbad_PNAS'],'epsc')
end
end