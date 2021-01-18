clear,clc,close all
folderName = '/home/crichter/ccsl/';
spreadsheet = dlmread('/home/crichter/ccsl/selecteddata_defcorr_all_70.csv',',',2,0);

% Capture only the withheld experiments for testing
exps = [8 26 33 37 39 48 52 53 95 103]; % 1,78 Removed
EXPERIMENT_USED_TO_SHOWCASE_EQ4_COMPONENTS = 16; %%%%%%%%%%%

totalexps = unique(spreadsheet(:,1));
exps = totalexps(totalexps>0);
nExps = max(exps);
eqdataExps = cell(nExps,1);
for i = 1:length(exps)
    whichExp = unique(spreadsheet(spreadsheet(:,1)==exps(i),1));
    eqdataExps{whichExp} = spreadsheet(spreadsheet(:,1)==whichExp,2:end);
end

toplot = [16 33 48 39 16];

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

thickness = .000096; % m
rho720 = 1215.9; % kg/m^3
##rho720 = 1215.92; % kg/m^3 <<-- Changed to 1215.9 to be consistent with other fn's.
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

%% RJW Calculation
Clmax = 1.9;
beta_trans = 1;
beta_rot = 1;

% Expressions for wy,wz taken from RJW equation 2.8, page 200
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

##Z01 = beta_trans*-lambdaz*.5.*(wydot-wx.*wz).*-(outerspan.^2-innerspan.^2);
##Z02 = beta_rot*-lambdazw.*wxdot.*(outerspan-innerspan);
##Z0spanwise = Z01+Z02;
##AMrjw = Z0spanwise.*-sin(def);
##RJWtot = Lrjw+AMrjw+inertial+inertialdev;

eq1 = 1.4526972e-8.*1000.*span.*Ycm.*defdot.*defdot.*cos(def);
C2 = 1.2489714e-7;
eq2 = C2.*1000.*strokedot.*span.*Xcm.*def;
eq3 = -1.0872654e-6.*1000.*span.*devddot;

beta_eq4 = [1.4526972e-8 1.2489714e-7 1.0872654e-6];
eq4_term_1 = fn_eq4_term_1(beta_eq4(1),eqdata);
eq4_term_2 = fn_eq4_term_2(beta_eq4(2),eqdata);
eq4_term_3 = fn_eq4_term_3(beta_eq4(3),eqdata);

##max(abs(eq4_term_1 - eq1))
##max(abs(eq4_term_2 - eq2))
##max(abs(eq4_term_3 - eq3))

##output = C1.*span.*Ycm.*defdot.*defdot.*cos(def)...
##    + C2.*strokedot.*span.*Xcm.*def...
##    - C3.*span.*devddot;

% Overfit equation
##eqbad = (14.9488 + 13.5087*def.*def - 1.38486.*def)./(407.258 + devddot + 17.8828.*devdot + cos(2.69241 - 42.9085./(13.5087.*def.*def) - 163.587.*def)) - 0.0381286;
beta_eqbad = [14.9488, 13.5087, 1.38486, 407.258, 17.8828, 2.69241, 42.9085, 13.5087, 163.587, 0.0381286];
eqbad = fn_eqbad(beta_eqbad, eqdata);

%----

%% Plotting Constants
domain = 0:length(force)-1;
axbuf = 1.2;
shadecolor = .9*[1 1 1];

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
           
figure(1)
set(gcf,'units','centimeters','paperunits','centimeters','papersize',[widthcm heightcm],...
    'paperposition',[0 0 widthcm heightcm])
set(gcf,'defaultlinelinewidth',.5,'visible','on')

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
saveFolder = '/home/crichter/ccsl/RSI21_Richter/paper/figures/';
##saveas(gcf,[saveFolder,'eq_comparison_PNAS'],'epsc')
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
set(gcf,'defaultlinelinewidth',.5,'visible','on')

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

if bla == 4
saveFolder = '/home/crichter/ccsl/RSI21_Richter/paper/figures/';
##saveas(gcf,[saveFolder,'eq_goodbad_PNAS'],'epsc')
end
end