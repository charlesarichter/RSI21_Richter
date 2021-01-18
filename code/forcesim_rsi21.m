clear,clc,close all
folderName = '/home/crichter/ccsl/';
spreadsheet = dlmread('/home/crichter/ccsl/selecteddata_defcorr_all_70.csv',',',2,0);
saveFolder = '/home/crichter/ccsl/RSI21_Richter/paper/figures/';

%% Compare terms of EQ4 with analytical model components
n_exp = 16;
eqdata = spreadsheet(spreadsheet(:,1)==n_exp,2:end);

%% Unpack Data
force = eqdata(:,1);
strokedot = eqdata(:,3);

%% RJW Calculation
Clmax = 1.9;
[~,Lrjw,inertial,inertialdev] = fn_rw_liftonly(Clmax,eqdata);

%% EQ4 (individual terms) Calculations
beta_eq4 = [1.4526972e-8 1.2489714e-7 1.0872654e-6];
[~,eq1,eq2,eq3] = fn_eq4(beta_eq4,eqdata);

%% Plotting Constants
axbuf = 1.2;
shadecolor = .9*[1 1 1];
fontsz = 12;
plotfactor = 1000; % Multiply force by 1000 then explain that all forces are "e-3"

cross1 = find([0;diff(sign(strokedot))]);
normstroke = [1;cross1;length(force)];
xticklab = {'0','.5','1.0','1.5','2'};

width = 15;
height = 7;

f1 = figure('visible','off');hold on
set(f1, 'defaultlinelinewidth',1)
fill([0 cross1(1) cross1(1) 0],10000*[-1 -1 1 1],shadecolor,'edgecolor',shadecolor)
fill([cross1(2) cross1(3) cross1(3) cross1(2)],10000*[-1 -1 1 1],shadecolor,'edgecolor',shadecolor)
plot(force*plotfactor,'k')
plot(eq2*plotfactor,'b')
plot(Lrjw*plotfactor,'r')

ylimit = axbuf*[min([force*plotfactor;(eq1+eq3)*plotfactor]) max([force*plotfactor;(eq1+eq3)*plotfactor])];
ylim(ylimit)
xlim([1 length(force)])
set(f1,'units','centimeters','paperunits','centimeters','papersize',[width height],...
    'paperposition',[0 0 width height])
set(get(f1,'currentaxes'),'box','on','layer','top','Xtick',normstroke,'XTickLabel',xticklab,'FontUnits','points','FontSize',fontsz)
xlabel('Stroke Cycle')
ylabel('Force (10^-^3 N)','FontSize',fontsz)
saveas(gcf,[saveFolder,'eq_comparison_lift.eps'],'epsc')

f2 = figure('visible','off');hold on
set(f1, 'defaultlinelinewidth',1)
fill([0 cross1(1) cross1(1) 0],10000*[-1 -1 1 1],shadecolor,'edgecolor',shadecolor)
fill([cross1(2) cross1(3) cross1(3) cross1(2)],10000*[-1 -1 1 1],shadecolor,'edgecolor',shadecolor)
plot(force*plotfactor,'k')
plot((eq1+eq3)*plotfactor,'b')
plot((inertial+inertialdev)*plotfactor,'r')

ylimit = axbuf*[min([force*plotfactor;(eq1+eq3)*plotfactor]) max([force*plotfactor;(eq1+eq3)*plotfactor])];
ylim(ylimit)
xlim([1 length(force)])
set(f2,'units','centimeters','paperunits','centimeters','papersize',[width height],...
    'paperposition',[0 0 width height])
set(get(f2,'currentaxes'),'box','on','layer','top','Xtick',normstroke,'XTickLabel',xticklab,'FontUnits','points','FontSize',fontsz)
xlabel('Stroke Cycle')
ylabel('Force (10^-^3 N)','FontSize',fontsz)
saveas(gcf,[saveFolder,'eq_comparison_inertial.eps'],'epsc')


%% Plot comparison between EQ4 and overfit Eureqa equation
toplot = [16 33 48 39];
for i_exp = 1:length(toplot)
eqdata = spreadsheet(spreadsheet(:,1)==toplot(i_exp),2:end);

%% Unpack Data
force = eqdata(:,1);
strokedot = eqdata(:,3);

%% EQ4 (individual terms) Calculations
beta_eq4 = [1.4526972e-8 1.2489714e-7 1.0872654e-6];
[~,eq1,eq2,eq3] = fn_eq4(beta_eq4,eqdata);

% Overfit equation
beta_eqbad = [14.9488, 13.5087, 1.38486, 407.258, 17.8828, 2.69241, 42.9085, 13.5087, 163.587, 0.0381286];
eqbad = fn_eqbad(beta_eqbad, eqdata);

%% Plotting Constants
domain = 0:length(force)-1;
axbuf = 1.2;
shadecolor = .9*[1 1 1];

fontsz = 12;
widthcm2 = 15;
heightcm2 = 9;
plotfactor = 1000; % Multiply force by 1000 then explain that all forces are "e-3"
xticklab = {'0','.5','1.0','1.5','2'};
        
f = figure('visible','off');hold on
set(f,'units','centimeters','paperunits','centimeters','papersize',[widthcm2 heightcm2],...
    'paperposition',[0 0 widthcm2 heightcm2])
set(f,'defaultlinelinewidth',1)

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

ylimit = axbuf*[min([force*plotfactor;(eqbad)*plotfactor;(eq1+eq2+eq3)*plotfactor])...
    max([force*plotfactor;(eqbad)*plotfactor;(eq1+eq2+eq3)*plotfactor])];
ylim(ylimit)
xlim([1 length(force)])
set(gca,'box','on','layer','top','Xtick',normstroke,'XTickLabel',xticklab,'FontUnits','points','FontSize',fontsz)
ylabel('Force (10^-^3 N)')
xlabel('Stroke Cycle')

saveas(f,[saveFolder,'eq_goodbad_',int2str(i_exp),'.eps'],'epsc')
end