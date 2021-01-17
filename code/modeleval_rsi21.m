clear,clc,close all
pkg load optim
pkg load signal
addpath("/home/crichter/ccsl/RSI21_Richter/code/models")

savefolder = '/home/crichter/ccsl/RSI21_Richter/paper/figures/';

%% Load Data
spreadsheet = dlmread('/home/crichter/ccsl/selecteddata_defcorr_all_70.csv',',',2,0);

%% Determine training, testing and other removed experiment IDs
exps_all = unique(spreadsheet(:,1));
exps_testing = [8 26 33 37 39 48 52 53 95 103]; % 1,78 removed
exps_removed = [1, 78];
exps_training = setdiff(exps_all, [exps_testing, exps_removed, 0]);

##% Mask off data to only include rows for which the experiment number (first column) is not in exps_testing, exps_removed, or 0.
##exps_training_mask = any(repmat(spreadsheet(:,1),1,size(exps_training)(1)) == repmat(exps_training',size(spreadsheet)(1),1),2);
##eqdata = spreadsheet(exps_training_mask,2:end);
##force = spreadsheet(exps_training_mask,2);

%% Equation sizes
sz_plpi = 32;
sz_plti = 46;
sz_rw_liftonly = 96;
sz_rw_liftam = 121;
sz_jw = 130;
sz_trans_rot = 118;
sz_trans_amrot = 113;
sz_trans_amtrans = 115;
sz_trans_amtrans_amrot = 120;
sz_trans_amrotddot = 106;
sz_eq1 = 24;
sz_eq2 = 19;
sz_eq3 = 16;
sz_eq4 = 14;
sz_eq5 = 13;
sz_eq6 = 12;

% Planar Inertial -- fn_planarinertial -- COMPLEXITY 13
% Total Inertial -- fn_totalinertial -- COMPLEXITY 27
% Planar Lift Planar Inertial -- fn_planarlift_planarinertial -- COMPLEXITY 32
% Planar Lift Total Inertial -- fn_planarlift_totalinertial -- COMPLEXITY 46
% RW Lift Only Total Inertial -- fn_rw_liftonly -- COMPLEXITY 96
% RW Lift and AM Total Inertial -- fn_rw_liftam -- COMPLEXITY 121
% JW Lift and Total Inertial -- fn_jw -- COMPLEXITY 130
% Analytical: fn_trans_rot -- fn_trans_rot -- COMPLEXITY 118
% Analytical: fn_trans_amrot -- fn_trans_amrot -- COMPLEXITY 113
% Analytical: fn_trans_amtrans -- fn_trans_amtrans -- COMPLEXITY 115
% Analytical: fn_trans_amtrans_amrot -- fn_trans_amtrans_amrot -- COMPLEXITY 120
% Analytical: fn_trans_amrotddot -- fn_trans_amrotddot -- COMPLEXITY 106
% EQ 1 -- fn_eq1 -- COMPLEXITY 24
% EQ 2 -- fn_eq2 -- COMPLEXITY 19
% EQ 3 -- fn_eq3 -- COMPLEXITY 16
% EQ 4 -- fn_eq4 -- COMPLEXITY 14
% EQ 5 -- fn_eq5 -- COMPLEXITY 13
% EQ 6 -- fn_eq6 -- COMPLEXITY 12

%% Coefficients of Eureqa models (fit by Eureqa to training data)
beta_eq1nf = [1.5723863e-8 1.2146853e-7 1.2146853e-7 1.2146853e-7 7.0700918e-7];
beta_eq2nf = [1.419276e-8 1.2280232e-7 1.2280232e-7 1.1480129e-6];
beta_eq3nf = [1.4180717e-8 1.2409225e-7 1.0660193e-6];
beta_eq4nf = [1.4526972e-8 1.2489714e-7 1.0872654e-6];
beta_eq5nf = [1.2524723e-7 5.8678208e-11 1.1539232e-6];
beta_eq6nf = [1.3552796e-8 1.2505514e-7 1.1357268e-6];

%% Specify initial values for parameter fitting
beta0_plpi = 1.5;
beta0_plti = 1.5;
beta0_rw_liftonly = 1.5;
beta0_rw_liftam = [1.5 1 1];
beta0_jw = [1.5 pi 1];
beta0_trans_rot = [1.5 2];
beta0_trans_amrot = [1.5 1];
beta0_trans_amtrans = [1.5 1];
beta0_trans_amtrans_amrot = [1.5 1 1];
beta0_trans_amrotddot = [1.5 1];

%% Fit coefficients of analyical models on all training data
disp("Fitting coefficients to training data...");
for i = 1:length(exps_training)
  i_exp = exps_training(i);
  eqdata = spreadsheet(spreadsheet(:,1)==i_exp,2:end);
  force = eqdata(:,1);
  
  beta_plpi(i,:) = nlinfit(eqdata,force,@fn_planarlift_planarinertial,beta0_plpi);
  beta_plti(i,:) = nlinfit(eqdata,force,@fn_planarlift_totalinertial,beta0_plti);
  beta_rw_liftonly(i,:) = nlinfit(eqdata,force,@fn_rw_liftonly,beta0_rw_liftonly);
  beta_rw_liftam(i,:) = nlinfit(eqdata,force,@fn_rw_liftam,beta0_rw_liftam);
  beta_jw(i,:) = nlinfit(eqdata,force,@fn_jw,beta0_jw);
  beta_trans_rot(i,:) = nlinfit(eqdata,force,@fn_trans_rot,beta0_trans_rot);
  beta_trans_amrot(i,:) = nlinfit(eqdata,force,@fn_trans_amrot,beta0_trans_amrot);
  beta_trans_amtrans(i,:) = nlinfit(eqdata,force,@fn_trans_amtrans,beta0_trans_amtrans);
  beta_trans_amtrans_amrot(i,:) = nlinfit(eqdata,force,@fn_trans_amtrans_amrot,beta0_trans_amtrans_amrot);
  beta_trans_amrotddot(i,:) = nlinfit(eqdata,force,@fn_trans_amrotddot,beta0_trans_amrotddot);
  
  beta_eq1(i,:) = nlinfit(eqdata,force,@fn_eq1,beta_eq1nf);
  beta_eq2(i,:) = nlinfit(eqdata,force,@fn_eq2,beta_eq2nf);
  beta_eq3(i,:) = nlinfit(eqdata,force,@fn_eq3,beta_eq3nf);
  beta_eq4(i,:) = nlinfit(eqdata,force,@fn_eq4,beta_eq4nf);
  beta_eq5(i,:) = nlinfit(eqdata,force,@fn_eq5,beta_eq5nf);
  beta_eq6(i,:) = nlinfit(eqdata,force,@fn_eq6,beta_eq6nf);
 
end

%% Mean values of coefficients for analytical models fit to training data
beta_plpi_mean = mean(beta_plpi, 1);
beta_plti_mean = mean(beta_plti, 1);
beta_rw_liftonly_mean = mean(beta_rw_liftonly, 1);
beta_rw_liftam_mean = mean(beta_rw_liftam, 1);
beta_jw_mean = mean(beta_jw, 1);
beta_trans_rot_mean = mean(beta_trans_rot, 1);
beta_trans_amrot_mean = mean(beta_trans_amrot, 1);
beta_trans_amtrans_mean = mean(beta_trans_amtrans, 1);
beta_trans_amtrans_amrot_mean = mean(beta_trans_amtrans_amrot, 1);
beta_trans_amrotddot_mean = mean(beta_trans_amrotddot, 1);

beta_eq1_mean = mean(beta_eq1, 1);
beta_eq2_mean = mean(beta_eq2, 1);
beta_eq3_mean = mean(beta_eq3, 1);
beta_eq4_mean = mean(beta_eq4, 1);
beta_eq5_mean = mean(beta_eq5, 1);
beta_eq6_mean = mean(beta_eq6, 1);

clear beta_plpi;
clear beta_plti;
clear beta_rw_liftonly;
clear beta_rw_liftam;
clear beta_jw;
clear beta_trans_rot;
clear beta_trans_amrot;
clear beta_trans_amtrans;
clear beta_trans_amtrans_amrot;
clear beta_trans_amrotddot;
clear beta_eq1
clear beta_eq2
clear beta_eq3
clear beta_eq4
clear beta_eq5
clear beta_eq6

%% Test all models on test data using coefficients fit to training data
disp("Testing fitted coefficients on testing data...");
for i = 1:length(exps_testing)
  i_exp = exps_testing(i);
  eqdata = spreadsheet(spreadsheet(:,1)==i_exp,2:end);
  force = eqdata(:,1);

  nofit_mae_plpi(i) = mean(abs(force-fn_planarlift_planarinertial(beta_plpi_mean,eqdata)));
  nofit_mae_plti(i) = mean(abs(force-fn_planarlift_totalinertial(beta_plti_mean,eqdata)));
  nofit_mae_rw_liftonly(i) = mean(abs(force-fn_rw_liftonly(beta_rw_liftonly_mean,eqdata)));
  nofit_mae_rw_liftam(i) = mean(abs(force-fn_rw_liftam(beta_rw_liftam_mean,eqdata)));
  nofit_mae_jw(i) = mean(abs(force-fn_jw(beta_jw_mean,eqdata)));
  nofit_mae_trans_rot(i) = mean(abs(force-fn_trans_rot(beta_trans_rot_mean,eqdata)));
  nofit_mae_trans_amrot(i) = mean(abs(force-fn_trans_amrot(beta_trans_amrot_mean,eqdata)));
  nofit_mae_trans_amtrans(i) = mean(abs(force-fn_trans_amtrans(beta_trans_amtrans_mean,eqdata)));
  nofit_mae_trans_amtrans_amrot(i) = mean(abs(force-fn_trans_amtrans_amrot(beta_trans_amtrans_amrot_mean,eqdata)));
  nofit_mae_trans_amrotddot(i) = mean(abs(force-fn_trans_amrotddot(beta_trans_amrotddot_mean,eqdata)));
  
  % Using coefficients fitted by Eureqa
  nofit_mae_eq1(i) = mean(abs(force-fn_eq1(beta_eq1nf,eqdata)));
  nofit_mae_eq2(i) = mean(abs(force-fn_eq2(beta_eq2nf,eqdata)));
  nofit_mae_eq3(i) = mean(abs(force-fn_eq3(beta_eq3nf,eqdata)));
  nofit_mae_eq4(i) = mean(abs(force-fn_eq4(beta_eq4nf,eqdata)));
  nofit_mae_eq5(i) = mean(abs(force-fn_eq5(beta_eq5nf,eqdata)));
  nofit_mae_eq6(i) = mean(abs(force-fn_eq6(beta_eq6nf,eqdata)));

  % Using coefficients fitted by nlinfit
##  mae_eq1(i) = mean(abs(force-fn_eq1(beta_eq1_mean,eqdata)));
##  mae_eq2(i) = mean(abs(force-fn_eq2(beta_eq2_mean,eqdata)));
##  mae_eq3(i) = mean(abs(force-fn_eq3(beta_eq3_mean,eqdata)));
##  mae_eq4(i) = mean(abs(force-fn_eq4(beta_eq4_mean,eqdata)));
##  mae_eq5(i) = mean(abs(force-fn_eq5(beta_eq5_mean,eqdata)));
##  mae_eq6(i) = mean(abs(force-fn_eq6(beta_eq6_mean,eqdata)));
end

format short e
disp("Analytical models for publication");
mean_nofit_mae_plpi = mean(nofit_mae_plpi,2)
mean_nofit_mae_plti = mean(nofit_mae_plti,2)
mean_nofit_mae_rw_liftonly = mean(nofit_mae_rw_liftonly,2)
mean_nofit_mae_rw_liftam = mean(nofit_mae_rw_liftam,2)
mean_nofit_mae_jw = mean(nofit_mae_jw,2)
disp("Other models");
mean_nofit_mae_trans_rot = mean(nofit_mae_trans_rot,2)
mean_nofit_mae_trans_amrot = mean(nofit_mae_trans_amrot,2)
mean_nofit_mae_trans_amtrans = mean(nofit_mae_trans_amtrans,2)
mean_nofit_mae_trans_amtrans_amrot = mean(nofit_mae_trans_amtrans_amrot,2)
mean_nofit_mae_trans_amrotddot = mean(nofit_mae_trans_amrotddot,2)
disp("Eureqa models");
mean_nofit_mae_eq1 = mean(nofit_mae_eq1,2)
mean_nofit_mae_eq2 = mean(nofit_mae_eq2,2)
mean_nofit_mae_eq3 = mean(nofit_mae_eq3,2)
mean_nofit_mae_eq4 = mean(nofit_mae_eq4,2)
##mean_nofit_mae_eq5 = mean(nofit_mae_eq5,2)
mean_nofit_mae_eq6 = mean(nofit_mae_eq6,2)
format

f = figure('visible','off');clf,hold on
largesize = 30;
smallsize = 20;
rotangle = 30;
fsize = 12;
##width=8.7;height=5;
width=16;height=10;
a = get(f,'currentaxes');
set(a,'box','on','layer','top','fontunits','points','fontsize',fsize)
set(f,'units','centimeters','paperunits','centimeters','papersize',[width height],...
    'paperposition',[0 0 width height])

plot(sz_plpi,mean_nofit_mae_plpi,'b.','Markersize',largesize);
text(sz_plpi,mean_nofit_mae_plpi,'   U^2 Planar 1','rotation',rotangle,'fontsize',fsize);
plot(sz_plti,mean_nofit_mae_plti,'g.','Markersize',largesize);
text(sz_plti,mean_nofit_mae_plti,'   U^2 Planar 2','rotation',rotangle,'fontsize',fsize);
plot(sz_rw_liftonly,mean_nofit_mae_rw_liftonly,'c.','Markersize',largesize);
text(sz_rw_liftonly,mean_nofit_mae_rw_liftonly,'   U^2 Total','rotation',rotangle,'fontsize',fsize);
plot(sz_rw_liftam,mean_nofit_mae_rw_liftam,'r.','Markersize',largesize);
text(sz_rw_liftam,mean_nofit_mae_rw_liftam,'   Whitney-Wood','rotation',rotangle,'fontsize',fsize);
plot(sz_jw,mean_nofit_mae_jw,'m.','Markersize',largesize);
text(sz_jw,mean_nofit_mae_jw,'   Pesavento-Wang','rotation',rotangle,'fontsize',fsize);
plot(sz_eq1,mean_nofit_mae_eq1,'k.','Markersize',smallsize);
plot(sz_eq2,mean_nofit_mae_eq2,'k.','Markersize',smallsize);
plot(sz_eq3,mean_nofit_mae_eq3,'k.','Markersize',smallsize);
plot(sz_eq4,mean_nofit_mae_eq4,'k.','Markersize',smallsize);
plot(sz_eq6,mean_nofit_mae_eq6,'k.','Markersize',smallsize);
text(sz_eq1-6,mean_nofit_mae_eq1-1.8e-4,{'EQ';'Models'},'horizontalalignment','center','fontsize',fsize);
ylabel('Mean Absolute Error (10^-^3 N)','fontsize',fsize)
xlabel('Equation Size (number of operators)','fontsize',fsize)

xmin = 0;
xmax = 180;
ymin = 0;
ymax = 1.9e-3;
xlim([xmin xmax])
ylim([ymin ymax])
dy = 2e-4;
yt = ymin:dy:ymax;
yticks(yt)
yticklabels(yt*1e3);
dx = 20;
xt = xmin:dx:xmax;
xticks(xt)
saveas(f,[savefolder, 'mae_nofit.eps'],'epsc');

%% Test all models on test data using coefficients fit to TEST data
disp("Fitting and testing on testing data...");
for i = 1:length(exps_testing)
  i_exp = exps_testing(i);
  eqdata = spreadsheet(spreadsheet(:,1)==i_exp,2:end);
  force = eqdata(:,1);

  beta_plpi = nlinfit(eqdata,force,@fn_planarlift_planarinertial,beta0_plpi);
  fit_mae_plpi(i) = mean(abs(force-fn_planarlift_planarinertial(beta_plpi,eqdata)));
  
  beta_plti = nlinfit(eqdata,force,@fn_planarlift_totalinertial,beta0_plti);
  fit_mae_plti(i) = mean(abs(force-fn_planarlift_totalinertial(beta_plti,eqdata)));
  
  beta_rw_liftonly = nlinfit(eqdata,force,@fn_rw_liftonly,beta0_rw_liftonly);
  fit_mae_rw_liftonly(i) = mean(abs(force-fn_rw_liftonly(beta_rw_liftonly,eqdata)));
  
  beta_rw_liftam = nlinfit(eqdata,force,@fn_rw_liftam,beta0_rw_liftam);
  fit_mae_rw_liftam(i) = mean(abs(force-fn_rw_liftam(beta_rw_liftam,eqdata)));
  
  beta_jw = nlinfit(eqdata,force,@fn_jw,beta0_jw);
  fit_mae_jw(i) = mean(abs(force-fn_jw(beta_jw,eqdata)));
  
  beta_trans_rot = nlinfit(eqdata,force,@fn_trans_rot,beta0_trans_rot);
  fit_mae_trans_rot(i) = mean(abs(force-fn_trans_rot(beta_trans_rot,eqdata)));

  beta_trans_amrot = nlinfit(eqdata,force,@fn_trans_amrot,beta0_trans_amrot);
  fit_mae_trans_amrot(i) = mean(abs(force-fn_trans_amrot(beta_trans_amrot,eqdata)));

  beta_trans_amtrans = nlinfit(eqdata,force,@fn_trans_amtrans,beta0_trans_amtrans);
  fit_mae_trans_amtrans(i) = mean(abs(force-fn_trans_amtrans(beta_trans_amtrans,eqdata)));

  beta_trans_amtrans_amrot = nlinfit(eqdata,force,@fn_trans_amtrans_amrot,beta0_trans_amtrans_amrot);
  fit_mae_trans_amtrans_amrot(i) = mean(abs(force-fn_trans_amtrans_amrot(beta_trans_amtrans_amrot,eqdata)));

  beta_trans_amrotddot = nlinfit(eqdata,force,@fn_trans_amrotddot,beta0_trans_amrotddot);
  fit_mae_trans_amrotddot(i) = mean(abs(force-fn_trans_amrotddot(beta_trans_amrotddot,eqdata)));
  
  beta_eq1 = nlinfit(eqdata,force,@fn_eq1,beta_eq1nf);
  fit_mae_eq1(i) = mean(abs(force-fn_eq1(beta_eq1,eqdata)));
  
  beta_eq2 = nlinfit(eqdata,force,@fn_eq2,beta_eq2nf);
  fit_mae_eq2(i) = mean(abs(force-fn_eq2(beta_eq2,eqdata)));
  
  beta_eq3 = nlinfit(eqdata,force,@fn_eq3,beta_eq3nf);
  fit_mae_eq3(i) = mean(abs(force-fn_eq3(beta_eq3,eqdata)));
  
  beta_eq4 = nlinfit(eqdata,force,@fn_eq4,beta_eq4nf);
  fit_mae_eq4(i) = mean(abs(force-fn_eq4(beta_eq4,eqdata)));
  
  beta_eq5 = nlinfit(eqdata,force,@fn_eq5,beta_eq5nf);
  fit_mae_eq5(i) = mean(abs(force-fn_eq5(beta_eq5,eqdata)));
  
  beta_eq6 = nlinfit(eqdata,force,@fn_eq6,beta_eq6nf);
  fit_mae_eq6(i) = mean(abs(force-fn_eq6(beta_eq6,eqdata)));
end

format short e
disp("Analytical models for publication");
mean(fit_mae_plpi,2)
mean(fit_mae_plti,2)
mean(fit_mae_rw_liftonly,2)
mean(fit_mae_rw_liftam,2)
mean(fit_mae_jw,2)
disp("Other models");
mean(fit_mae_trans_rot,2)
mean(fit_mae_trans_amrot,2)
mean(fit_mae_trans_amtrans,2)
mean(fit_mae_trans_amtrans_amrot,2)
mean(fit_mae_trans_amrotddot,2)
disp("Eureqa models");
mean(fit_mae_eq1,2)
mean(fit_mae_eq2,2)
mean(fit_mae_eq3,2)
mean(fit_mae_eq4,2)
##mean(fit_mae_eq5,2)
mean(fit_mae_eq6,2)
format

return


##nExps = length(exps);
##
##forcepredict = cell(1,nExps);
##forcepredicteq = cell(1,nExps);
##mae1 = nan(1,nExps);mae1eq = nan(1,nExps);
##rmd1 = nan(1,nExps);rmd1eq = nan(1,nExps);
##eqdata = cell(nExps,1);
##force = cell(nExps,1);
##
##% Initialize Coefficients - MEAN VALUES CALCULATED FROM *ALL* EXPS
##beta_jw = nan(nExps,3); % [Cl Crotlift Cm22] = [1.6864 1.4752 0.7536], STD = [0.2389 0.9920 0.4233]
##beta_rw_liftam = nan(nExps,3); % [Cl Camtrans Camrotddot]  = [1.8429 0.4382 -0.3152], STD = [0.2430 0.5158 1.0283]
##beta_trans_rot = nan(nExps,2); % [Cl Crotlift] = [1.7361 -0.0908], STD = [0.2453 0.8649]
##beta_trans_amrot = nan(nExps,2); % [Cl Camrot] = [1.7420 0.0755], STD = [0.2432 0.6126]
##beta_trans_amtrans = nan(nExps,2); % [Cl Camtrans] = [1.7354 0.3756], STD = [0.2514 0.5499]
##beta_trans_amtrans_amrot = nan(nExps,3); % [Cl Camtrans Camrot] = [1.7109 0.6684 -0.3172], STD = [0.2316 0.4530 0.5491]
##beta_trans_amtransRW = nan(nExps,2);
##beta_trans_amrotddot = nan(nExps,2); % [Cl Camrotddot] = [1.8360 -0.2577], STD = [0.2404 0.9987]
##
##% EQ coefficients
##beta_eq1nf = [1.5723863e-8 1.2146853e-7 1.2146853e-7 1.2146853e-7 7.0700918e-7];
##beta_eq2nf = [1.419276e-8 1.2280232e-7 1.2280232e-7 1.1480129e-6];
##beta_eq3nf = [1.4180717e-8 1.2409225e-7 1.0660193e-6];
##beta_eq4nf = [1.4526972e-8 1.2489714e-7 1.0872654e-6];
##beta_eq5nf = [1.2524723e-7 5.8678208e-11 1.1539232e-6];
##beta_eq6nf = [1.3552796e-8 1.2505514e-7 1.1357268e-6];
##
##%% Loop Through Test Experiments - With Coefficient Fitting
##for i = 1:nExps
##    disp(i)
##    eqdata{i} = spreadsheet(spreadsheet(:,1)==exps(i),2:end);
##    force{i} = eqdata{i}(:,1);
##
##    % Planar Inertial
##    forcepredict{1,i} = fn_planarinertial(eqdata{i});
##    mae1(1,i) = mean(abs(force{i}-forcepredict{1,i}));rmd1(1,i) = mae1(1,i)/mean(force{i});
##
##    % Total Inertial
##    forcepredict{2,i} = fn_totalinertial(eqdata{i});
##    mae1(2,i) = mean(abs(force{i}-forcepredict{2,i}));rmd1(2,i) = mae1(2,i)/mean(force{i});
##
##    % Planar Lift Planar Inertial
##    beta0_plpi = 1.5;
##    beta_plpi(i) = nlinfit(eqdata{i},force{i},@fn_planarlift_planarinertial,beta0_plpi); 
##    forcepredict{3,i} = fn_planarlift_planarinertial(beta_plpi(i),eqdata{i});
##    mae1(3,i) = mean(abs(force{i}-forcepredict{3,i}));rmd1(3,i) = mae1(3,i)/mean(force{i});
##
##    % Planar Lift Total Inertial
##    beta0_plti = 1.5;
##    beta_plti(i) = nlinfit(eqdata{i},force{i},@fn_planarlift_totalinertial,beta0_plti); 
##    forcepredict{4,i} = fn_planarlift_totalinertial(beta_plti(i),eqdata{i});
##    mae1(4,i) = mean(abs(force{i}-forcepredict{4,i}));rmd1(4,i) = mae1(4,i)/mean(force{i});
##
##    % RW Lift Only Total Inertial
##    beta0_rw_liftonly = 1.5;
##    beta_rw_liftonly(i) = nlinfit(eqdata{i},force{i},@fn_rw_liftonly,beta0_rw_liftonly);
##    forcepredict{5,i} = fn_rw_liftonly(beta_rw_liftonly(i),eqdata{i});
##    mae1(5,i) = mean(abs(force{i}-forcepredict{5,i}));rmd1(5,i) = mae1(5,i)/mean(force{i});
##
##    % RW Lift and AM Total Inertial
##    beta0_rw_liftam = [1.5 1 1];
##    beta_rw_liftam(i,:) = nlinfit(eqdata{i},force{i},@fn_rw_liftam,beta0_rw_liftam);
##    forcepredict{6,i} = fn_rw_liftam(beta_rw_liftam(i,:),eqdata{i});
##    mae1(6,i) = mean(abs(force{i}-forcepredict{6,i}));rmd1(6,i) = mae1(6,i)/mean(force{i});
##
##    % JW Lift and Total Inertial
##    beta0_jw = [1.5 pi 1];
##    beta_jw(i,:) = nlinfit(eqdata{i},force{i},@fn_jw,beta0_jw);
##    forcepredict{7,i} = fn_jw(beta_jw(i,:),eqdata{i});
##    mae1(7,i) = mean(abs(force{i}-forcepredict{7,i}));rmd1(7,i) = mae1(7,i)/mean(force{i});
##
##    % Analytical: fn_trans_rot -- COMPLEXITY 118
##    beta0_trans_rot = [1.5 2];
##    beta_trans_rot(i,:) = nlinfit(eqdata{i},force{i},@fn_trans_rot,beta0_trans_rot);
##    forcepredict{8,i} = fn_trans_rot(beta_trans_rot(i,:),eqdata{i});
##    mae1(8,i) = mean(abs(force{i}-forcepredict{8,i}));rmd1(8,i) = mae1(8,i)/mean(force{i});
##    corrholder = corrcoef(force{i},forcepredict{8,i});corr1(8,i) = corrholder(2);
##
##    % Analytical: fn_trans_amrot -- COMPLEXITY 113
##    beta0_trans_amrot = [1.5 1];
##    beta_trans_amrot(i,:) = nlinfit(eqdata{i},force{i},@fn_trans_amrot,beta0_trans_amrot);
##    forcepredict{9,i} = fn_trans_amrot(beta_trans_amrot(i,:),eqdata{i});
##    mae1(9,i) = mean(abs(force{i}-forcepredict{9,i}));rmd1(9,i) = mae1(9,i)/mean(force{i});
##    corrholder = corrcoef(force{i},forcepredict{9,i});corr1(9,i) = corrholder(2);
##
##    % Analytical: fn_trans_amtrans -- COMPLEXITY 115
##    beta0_trans_amtrans = [1.5 1];
##    beta_trans_amtrans(i,:) = nlinfit(eqdata{i},force{i},@fn_trans_amtrans,beta0_trans_amtrans);
##    forcepredict{10,i} = fn_trans_amtrans(beta_trans_amtrans(i,:),eqdata{i});
##    mae1(10,i) = mean(abs(force{i}-forcepredict{10,i}));rmd1(10,i) = mae1(10,i)/mean(force{i});
##    corrholder = corrcoef(force{i},forcepredict{10,i});corr1(10,i) = corrholder(2);
##
##    % Analytical: fn_trans_amtrans_amrot -- COMPLEXITY 120
##    beta0_trans_amtrans_amrot = [1.5 1 1];
##    beta_trans_amtrans_amrot(i,:) = nlinfit(eqdata{i},force{i},@fn_trans_amtrans_amrot,beta0_trans_amtrans_amrot);
##    forcepredict{11,i} = fn_trans_amtrans_amrot(beta_trans_amtrans_amrot(i,:),eqdata{i});
##    mae1(11,i) = mean(abs(force{i}-forcepredict{11,i}));rmd1(11,i) = mae1(11,i)/mean(force{i});
##    corrholder = corrcoef(force{i},forcepredict{11,i});corr1(11,i) = corrholder(2);
##
##    % Analytical: fn_trans_amrotddot -- COMPLEXITY 106
##    beta0_trans_amrotddot = [1.5 1];
##    beta_trans_amrotddot(i,:) = nlinfit(eqdata{i},force{i},@fn_trans_amrotddot,beta0_trans_amrotddot);
##    forcepredict{12,i} = fn_trans_amrotddot(beta_trans_amrotddot(i,:),eqdata{i});
##    mae1(12,i) = mean(abs(force{i}-forcepredict{12,i}));rmd1(12,i) = mae1(12,i)/mean(force{i});
##    corrholder = corrcoef(force{i},forcepredict{12,i});corr1(12,i) = corrholder(2);
##
##    % EQ 1
##    beta0_eq1 = [1e-7 1e-7 1e-7 1e-7 1e-7];
##    beta_eq1(i,:) = nlinfit(eqdata{i},force{i},@fn_eq1,beta0_eq1);
##    forcepredicteq{1,i} = fn_eq1(beta_eq1(i,:),eqdata{i});
##    mae1eq(1,i) = mean(abs(force{i}-forcepredicteq{1,i}));rmd1eq(1,i) = mae1eq(1,i)/mean(force{i});
##
##    % EQ 2
##    beta0_eq2 = [1e-7 1e-7 1e-7 1e-7];
##    beta_eq2(i,:) = nlinfit(eqdata{i},force{i},@fn_eq2,beta0_eq2);
##    forcepredicteq{2,i} = fn_eq2(beta_eq2(i,:),eqdata{i});
##    mae1eq(2,i) = mean(abs(force{i}-forcepredicteq{2,i}));rmd1eq(2,i) = mae1eq(2,i)/mean(force{i});
##
##    % EQ 3
##    beta0_eq3 = [1e-7 1e-7 1e-7];
##    beta_eq3(i,:) = nlinfit(eqdata{i},force{i},@fn_eq3,beta0_eq3);
##    forcepredicteq{3,i} = fn_eq3(beta_eq3(i,:),eqdata{i});
##    mae1eq(3,i) = mean(abs(force{i}-forcepredicteq{3,i}));rmd1eq(3,i) = mae1eq(3,i)/mean(force{i});
##    
##    % EQ 4
##    beta0_eq4 = [1e-7 1e-7 1e-7];
##    beta_eq4(i,:) = nlinfit(eqdata{i},force{i},@fn_eq4,beta0_eq4);
##    forcepredicteq{4,i} = fn_eq4(beta_eq4(i,:),eqdata{i});
##    mae1eq(4,i) = mean(abs(force{i}-forcepredicteq{4,i}));rmd1eq(4,i) = mae1eq(4,i)/mean(force{i});
##    
##    % EQ 5
##    beta0_eq5 = [1e-7 1e-7 1e-7];
##    beta_eq5(i,:) = nlinfit(eqdata{i},force{i},@fn_eq5,beta0_eq5);
##    forcepredicteq{5,i} = fn_eq5(beta_eq5(i,:),eqdata{i});
##    mae1eq(5,i) = mean(abs(force{i}-forcepredicteq{5,i}));rmd1eq(5,i) = mae1eq(5,i)/mean(force{i});
##    
##    % EQ 6
##    beta0_eq6 = [1e-7 1e-7 1e-7];
##    beta_eq6(i,:) = nlinfit(eqdata{i},force{i},@fn_eq6,beta0_eq6);
##    forcepredicteq{6,i} = fn_eq6(beta_eq6(i,:),eqdata{i});
##    mae1eq(6,i) = mean(abs(force{i}-forcepredicteq{6,i}));rmd1eq(6,i) = mae1eq(6,i)/mean(force{i});
##end
##%% Averages
##fitmae_an = mean(mae1,2);
##fitrmd_an = mean(rmd1,2);
##fitmae_eq = mean(mae1eq,2);
##fitrmd_eq = mean(rmd1eq,2);
##
##%% Coeff. Averages
##beta_plpi_nf = mean(beta_plpi);
##beta_plti_nf = mean(beta_plti);
##beta_rw_liftonly_nf = mean(beta_rw_liftonly);
##
##beta_rw_liftam_nf = mean(beta_rw_liftam,1);
##beta_rw_std = std(beta_rw_liftam,1);
##beta_jw_nf = mean(beta_jw,1);
##beta_jw_std = std(beta_jw,1);
##
##beta_trans_rot_nf = mean(beta_trans_rot,1);
##beta_trans_amrot_nf = mean(beta_trans_amrot,1);
##beta_trans_amtrans_nf = mean(beta_trans_amtrans,1);
##beta_trans_amtrans_amrot_nf = mean(beta_trans_amtrans_amrot,1);
##beta_trans_amrotddot_nf = mean(beta_trans_amrotddot,1);
##
##%% Coeff. Averages Jointly Optimized Over All Training Experiments
####beta_plpi_nf =  1.7154;
####beta_plti_nf =  1.6464;
####beta_rw_liftonly_nf =  1.6483;
####beta_rw_liftam_nf = [1.82200,  -0.13409,  -0.83757];
####beta_jw_nf =   [1.57661,   1.42240,   0.45386];
####beta_trans_rot_nf =   [1.61200,   0.63702];
####beta_trans_amrot_nf =   [1.62624,  -0.43529];
####beta_trans_amtrans_nf =   [1.65580,  -0.24950];
####beta_trans_amtrans_amrot_nf =   [1.61450,   0.21118,  -0.54186];
####beta_trans_amrotddot_nf =   [1.82460,  -0.87030];
##
##clear mae1 mae1eq rmd1 rmd1eq eqdata force
##
##%% Loop Through Test Experiments - No Coefficient Fitting
##forcepredict2 = cell(1,nExps);
##forcepredicteq2 = cell(1,nExps);
##mae2 = nan(1,nExps);mae2eq = nan(1,nExps);
##rmd2 = nan(1,nExps);rmd2eq = nan(1,nExps);
##eqdata = cell(nExps,1);
##force = cell(nExps,1);
##
##for i = 1:nExps
##    eqdata{i} = spreadsheet(spreadsheet(:,1)==exps(i),2:end);
##    force{i} = eqdata{i}(:,1);
##
##    % Planar Inertial -- COMPLEXITY 13
##    forcepredict2{1,i} = fn_planarinertial(eqdata{i});
##    mae2(1,i) = mean(abs(force{i}-forcepredict2{1,i}));rmd2(1,i) = mae2(1,i)/mean(force{i});
##    
##    % Total Inertial -- COMPLEXITY 27
##    forcepredict2{2,i} = fn_totalinertial(eqdata{i});
##    mae2(2,i) = mean(abs(force{i}-forcepredict2{2,i}));rmd2(2,i) = mae2(2,i)/mean(force{i});
##
##    % Planar Lift Planar Inertial -- COMPLEXITY 32
##%     beta_plpi_nf = 1.7942;
##    forcepredict2{3,i} = fn_planarlift_planarinertial(beta_plpi_nf,eqdata{i});
##    mae2(3,i) = mean(abs(force{i}-forcepredict2{3,i}));rmd2(3,i) = mae2(3,i)/mean(force{i});
##
##    % Planar Lift Total Inertial -- COMPLEXITY 46
##%     beta_plti_nf = 1.7327;
##    forcepredict2{4,i} = fn_planarlift_totalinertial(beta_plti_nf,eqdata{i});
##    mae2(4,i) = mean(abs(force{i}-forcepredict2{4,i}));rmd2(4,i) = mae2(4,i)/mean(force{i});
##
##    % RW Lift Only Total Inertial -- COMPLEXITY 96
##%     beta_rw_liftonly_nf = 1.7410;
##    forcepredict2{5,i} = fn_rw_liftonly(beta_rw_liftonly_nf,eqdata{i});
##    mae2(5,i) = mean(abs(force{i}-forcepredict2{5,i}));rmd2(5,i) = mae2(5,i)/mean(force{i});
##
##    % RW Lift and AM Total Inertial -- COMPLEXITY 121
##%     beta_rw_liftam_nf = [1.8429 0.4382 -0.3152];
##    forcepredict2{6,i} = fn_rw_liftam(beta_rw_liftam_nf,eqdata{i});
##    mae2(6,i) = mean(abs(force{i}-forcepredict2{6,i}));rmd2(6,i) = mae2(6,i)/mean(force{i});
##
##    % JW Lift and Total Inertial -- COMPLEXITY 130
##%     beta_jw_nf = [1.6864 1.4752 0.7536];
##    forcepredict2{7,i} = fn_jw(beta_jw_nf,eqdata{i});
##    mae2(7,i) = mean(abs(force{i}-forcepredict2{7,i}));rmd2(7,i) = mae2(7,i)/mean(force{i});
##
##    % Analytical: fn_trans_rot -- COMPLEXITY 118
##%     beta_trans_rot_nf = [1.7361 -0.0908];
##    forcepredict2{8,i} = fn_trans_rot(beta_trans_rot_nf,eqdata{i});
##    mae2(8,i) = mean(abs(force{i}-forcepredict2{8,i}));rmd2(8,i) = mae2(8,i)/mean(force{i});
##
##    % Analytical: fn_trans_amrot -- COMPLEXITY 113
##%     beta_trans_amrot_nf = [1.7361 -0.0908];
##    forcepredict2{9,i} = fn_trans_amrot(beta_trans_amrot_nf,eqdata{i});
##    mae2(9,i) = mean(abs(force{i}-forcepredict2{9,i}));rmd2(9,i) = mae2(9,i)/mean(force{i});
##
##    % Analytical: fn_trans_amtrans -- COMPLEXITY 115
##%     beta_trans_amtrans_nf = [1.7420 0.0755];
##    forcepredict2{10,i} = fn_trans_amtrans(beta_trans_amtrans_nf,eqdata{i});
##    mae2(10,i) = mean(abs(force{i}-forcepredict2{10,i}));rmd2(10,i) = mae2(10,i)/mean(force{i});
##
##    % Analytical: fn_trans_amtrans_amrot -- COMPLEXITY 120
##%     beta_trans_amtrans_amrot_nf = [1.7109 0.6684 -0.3172];
##    forcepredict2{11,i} = fn_trans_amtrans_amrot(beta_trans_amtrans_amrot_nf,eqdata{i});
##    mae2(11,i) = mean(abs(force{i}-forcepredict2{11,i}));rmd2(11,i) = mae2(11,i)/mean(force{i});
##
##    % Analytical: fn_trans_amrotddot -- COMPLEXITY 106
##%     beta_trans_amrotddot_nf = [1.8360 -0.2577];
##    forcepredict2{12,i} = fn_trans_amrotddot(beta_trans_amrotddot_nf,eqdata{i});
##    mae2(12,i) = mean(abs(force{i}-forcepredict2{12,i}));rmd2(12,i) = mae2(12,i)/mean(force{i});
##
##    % EQ 1 -- COMPLEXITY 24
##    forcepredicteq2{1,i} = fn_eq1(beta_eq1nf,eqdata{i});
##    mae2eq(1,i) = mean(abs(force{i}-forcepredicteq2{1,i}));rmd2eq(1,i) = mae2eq(1,i)/mean(force{i});
##    
##    % EQ 2 -- COMPLEXITY 19
##    forcepredicteq2{2,i} = fn_eq2(beta_eq2nf,eqdata{i});
##    mae2eq(2,i) = mean(abs(force{i}-forcepredicteq2{2,i}));rmd2eq(2,i) = mae2eq(2,i)/mean(force{i});
##    
##    % EQ 3 -- COMPLEXITY 16
##    forcepredicteq2{3,i} = fn_eq3(beta_eq3nf,eqdata{i});
##    mae2eq(3,i) = mean(abs(force{i}-forcepredicteq2{3,i}));rmd2eq(3,i) = mae2eq(3,i)/mean(force{i});
##    
##    % EQ 4 -- COMPLEXITY 14
##    forcepredicteq2{4,i} = fn_eq4(beta_eq4nf,eqdata{i});
##    mae2eq(4,i) = mean(abs(force{i}-forcepredicteq2{4,i}));rmd2eq(4,i) = mae2eq(4,i)/mean(force{i});
##    
##    % EQ 5 -- COMPLEXITY 13
##    forcepredicteq2{5,i} = fn_eq5(beta_eq5nf,eqdata{i});
##    mae2eq(5,i) = mean(abs(force{i}-forcepredicteq2{5,i}));rmd2eq(5,i) = mae2eq(5,i)/mean(force{i});
##    
##    % EQ 6 -- COMPLEXITY 12
##    forcepredicteq2{6,i} = fn_eq6(beta_eq6nf,eqdata{i});
##    mae2eq(6,i) = mean(abs(force{i}-forcepredicteq2{6,i}));rmd2eq(6,i) = mae2eq(6,i)/mean(force{i});
##end
##%% Averages
##nofitmae_an = mean(mae2,2);
##nofitrmd_an = mean(rmd2,2);
##nofitmae_eq = mean(mae2eq,2);
##nofitrmd_eq = mean(rmd2eq,2);
##
##%% Polynomial
##pOrder = 2:1:30;
##for n = 1:length(pOrder)
##    for k = 1:nExps
##        [pks,locs]=findpeaks(abs(eqdata{k}(:,2)),'minpeakheight',0);
##        middle = round(length(eqdata{k}(:,2))/2);        
##        [middleval,middlepos] = min(abs(locs-middle));
##        fpoly{k} = force{k}(locs(middlepos-1):locs(middlepos+1));
##        
##        polyx{k} = (1:length(fpoly{k}))';
##        [p{k},s{k},mu{k}] = polyfit(polyx{k},fpoly{k},pOrder(n));
##        y{k} = polyval(p{k},polyx{k},[],mu{k});
##        e{k} = fpoly{k}-y{k};
##        maepoly{k} = mean(abs(e{k}));
##    end
##    polymae(n) = mean(cell2mat(maepoly));
##
##%     figure(4)
##%     for m = 1:nExps
##%         subplot(nExps,1,m),hold on
##%         plot(polyx{m},fpoly{m},'k','linewidth',2)
##%         plot(polyx{m},y{m},'b')
##% %         plot(polyx{m},e{m},'g')
##%         ylim([-.015 .015]),hold off
##%     end
##%     figure(3),hold on
##%     plot(3*n-1,polymae(n),'k+')
##end
##% plot(3*pOrder-1,polymae,'k','linewidth',2)
##
##%% Complexity
##complex_an = [13 27 32 46 96 121 130 118 113 115 120 106];
##complex_eq = [24 19 16 14 13 12];
##
##%% Pareto - Fitted Coefficients
##figure(1),clf,hold on
##% set(gcf,'units','inches','outerposition',[0 0 5.8 3])
##smallmsize = 16;
##rotangle = 30;
##fsize = 6;
##
##%==== FOR PNAS - SWITCHED FROM NOFITMAE TO FITMAE
##% plot(complex_an(3),fitmae_an(3),'b.','Markersize',20)%smallmsize)
##% text(complex_an(3),fitmae_an(3),'   U^2 Planar 1','rotation',rotangle,'fontsize',fsize)
##% plot(complex_an(4),fitmae_an(4),'g.','Markersize',20)%,smallmsize)
##% text(complex_an(4),fitmae_an(4),'   U^2 Planar 2','rotation',rotangle,'fontsize',fsize)
##% 
##% plot(complex_an(5),fitmae_an(5),'c.','Markersize',20)%,smallmsize)
##% text(complex_an(5),fitmae_an(5),'   U^2 Total','rotation',rotangle,'fontsize',fsize)
##% 
##% % plot(complex_an(8:12),fitmae_an(8:12),'g.','Markersize',smallmsize)
##% % text(complex_an(11)-4,fitmae_an(11)+4e-5,'   Other Combinations','rotation',rotangle,'fontsize',fsize)
##% 
##% plot(complex_an(6),fitmae_an(6),'r.','Markersize',20)
##% text(complex_an(6),fitmae_an(6)+2e-5,'   Whitney-Wood','rotation',rotangle,'fontsize',fsize)
##% 
##% plot(complex_an(7),fitmae_an(7),'m.','Markersize',20)
##% text(complex_an(7),fitmae_an(7),'   Pesavento-Wang','rotation',rotangle,'fontsize',fsize)
##% 
##% plot(complex_eq([1,2,3,4,6]),fitmae_eq([1,2,3,4,6]),'k.','Markersize',smallmsize)
##% text(complex_eq(1)-6,fitmae_eq(1)-1.5e-4,{'EQ';'Models'},'horizontalalignment','center','fontsize',fsize)
##%====
##
##%==== 
##plot(complex_an(3),nofitmae_an(3),'b.','Markersize',20)%smallmsize)
##text(complex_an(3),nofitmae_an(3),'   U^2 Planar 1','rotation',rotangle,'fontsize',fsize)
##plot(complex_an(4),nofitmae_an(4),'g.','Markersize',20)%,smallmsize)
##text(complex_an(4),nofitmae_an(4),'   U^2 Planar 2','rotation',rotangle,'fontsize',fsize)
##
##plot(complex_an(5),nofitmae_an(5),'c.','Markersize',20)%,smallmsize)
##text(complex_an(5),nofitmae_an(5),'   U^2 Total','rotation',rotangle,'fontsize',fsize)
##
##% plot(complex_an(8:12),nofitmae_an(8:12),'g.','Markersize',smallmsize)
##% text(complex_an(11)-4,nofitmae_an(11)+4e-5,'   Other Combinations','rotation',rotangle,'fontsize',fsize)
##
##plot(complex_an(6),nofitmae_an(6),'r.','Markersize',20)
##text(complex_an(6),nofitmae_an(6)+2e-5,'   Whitney-Wood','rotation',rotangle,'fontsize',fsize)
##
##plot(complex_an(7),nofitmae_an(7),'m.','Markersize',20)
##text(complex_an(7),nofitmae_an(7),'   Pesavento-Wang','rotation',rotangle,'fontsize',fsize)
##
##plot(complex_eq([1,2,3,4,6]),nofitmae_eq([1,2,3,4,6]),'k.','Markersize',smallmsize)
##text(complex_eq(1)-6,nofitmae_eq(1)-1.5e-4,{'EQ';'Models'},'horizontalalignment','center','fontsize',fsize)
##
##% plot(complex_an(3:end),fitmae_an(3:end),'ro')
##% plot(complex_eq,fitmae_eq,'bo')
##% % % plot(3*pOrder-1,polymae,'k.','markersize',6)
##%====
##
##ylim([0 19e-4])
##xlim([0 180])
##ylabel('Mean Absolute Error (10^-^3 N)','fontsize',fsize)
##xlabel('Equation Size (number of operators)','fontsize',fsize)
##width=8.7;height=5;
##lbuf=1;rbuf=.2;
##bbuf=.8;tbuf=0;
####set(gca,'box','on','layer','top','fontunits','points','fontsize',fsize)
####set(gca,'units','centimeters','position',[lbuf bbuf width-lbuf-rbuf height-bbuf-tbuf])
####set(gcf,'units','centimeters','paperunits','centimeters','papersize',[width height],...
####    'paperposition',[0 0 width height])
##
##% saveFolder = '/Users/charlesrichter/Desktop/CR Shared Drive/Thesis/figures/';
####saveFolder = '/Users/charlesrichter/Desktop/PNAS11_Richter/figures/';
####saveFolder = '/home/crichter/ccsl/NSR20_Richter/paper/figures/';
##saveFolder = '/home/crichter/ccsl/RSI21_Richter/paper/figures/';
##saveas(gcf,[saveFolder,'pareto_mae_nofit_PNAS'],'epsc')
##% MARKERSIZES, FONTSIZES, DIMENSIONS, ETC ARE CHANGED FROM THE ORIGINAL
##% THESIS PLOT TO BE USED IN THE JOURNAL PUBLICATION.  REFER TO THE PLOTS
##% BELOW TO RETURN TO THE ORIGINAL VALUES.
##
##
##%% Pareto - Fitted Coefficients
##figure(2),clf,hold on
##% set(gcf,'units','inches','outerposition',[0 0 5.8 3])
##smallmsize = 16;
##rotangle = 30;
##fsize = 8;
##
##plot(complex_an(3),fitmae_an(3),'b.','Markersize',smallmsize)
##text(complex_an(3),fitmae_an(3),'  U^2 Planar 1','rotation',rotangle,'fontsize',fsize)
##plot(complex_an(4),fitmae_an(4),'b.','Markersize',smallmsize)
##text(complex_an(4),fitmae_an(4),'  U^2 Planar 2','rotation',rotangle,'fontsize',fsize)
##
##plot(complex_an(5),fitmae_an(5),'c.','Markersize',smallmsize)
##text(complex_an(5),fitmae_an(5),'  U^2 Total','rotation',rotangle,'fontsize',fsize)
##
##plot(complex_an(8:12),fitmae_an(8:12),'g.','Markersize',smallmsize)
##text(complex_an(10),fitmae_an(10)+2e-5,' Other Combinations','rotation',rotangle,'fontsize',fsize)
##
##plot(complex_an(6),fitmae_an(6),'r.','Markersize',24)
##text(complex_an(6),fitmae_an(6),'   Whitney-Wood','rotation',rotangle,'fontsize',fsize)
##
##plot(complex_an(7),fitmae_an(7),'m.','Markersize',24)
##text(complex_an(7),fitmae_an(7),'   Pesavento-Wang','rotation',rotangle,'fontsize',fsize)
##
##plot(complex_eq,fitmae_eq,'k.','Markersize',smallmsize)
##text(complex_eq(1)-6,fitmae_eq(1)-1.2e-4,{'EQ';'Models'},'horizontalalignment','center','fontsize',fsize)
##
##% plot(complex_an(3:end),fitmae_an(3:end),'ro')
##% plot(complex_eq,fitmae_eq,'bo')
##plot(3*pOrder-1,polymae,'k.','markersize',6)
##
##ylim([0 17e-4])
##xlim([0 170])
##ylabel('Mean Absolute Error')
##xlabel('Equation Size (number of operators)')
##width=5.8;height=3;
##    set(gca,'box','on','layer','top')
##set(gcf,'units','inches','papersize',[width height],'position',[1 9 width height],'paperpositionmode','auto')
##
##% saveFolder = '/Users/charlesrichter/Desktop/CR Shared Drive/Thesis/figures/';
##% saveas(gcf,[saveFolder,'pareto_mae_oct23'],'pdf')
##
##
##%% Force Plots
##figure(3),clf
##set(gcf,'defaultlinelinewidth',.5)
##iExp = [2 4 6 9];
##shadecolor = .9*[1,1,1];
##
##for n = 1:length(iExp)
##    subplot(2,2,n),hold on  
##    strokedot = eqdata{iExp(n)}(:,3);
##    cross1 = find([0;diff(sign(strokedot))]);
##    if cross1(1)<100
##        fill([cross1(1) cross1(2) cross1(2) cross1(1)],100*[-1 -1 1 1],shadecolor,'edgecolor',shadecolor)
##        fill([cross1(3) cross1(4) cross1(4) cross1(3)],100*[-1 -1 1 1],shadecolor,'edgecolor',shadecolor)
##        xlim([cross1(1) length(force{iExp(n)})])
##    else
##        fill([0 cross1(1) cross1(1) 0],100*[-1 -1 1 1],shadecolor,'edgecolor',shadecolor)
##        fill([cross1(2) cross1(3) cross1(3) cross1(2)],100*[-1 -1 1 1],shadecolor,'edgecolor',shadecolor)
##        xlim([1 length(force{iExp(n)})])
##    end
##    clear strokedot
##    plot(force{iExp(n)},'k')
##    ylimit = 1.2*[min(force{iExp(n)}) max(force{iExp(n)})];
##%     ylimit = 1.2*[min(force{iExp(n)}) max(forcepredict{6,iExp(n)})];
##
##    plot(forcepredict{3,iExp(n)},'b')
##    plot(forcepredict{4,iExp(n)},'r')
##
##%     plot(forcepredict{5,iExp(n)},'c')
##%     plot(forcepredict{6,iExp(n)},'m')
##
##%     plot(forcepredict{7,iExp(n)},'r')
##
##%     plot(forcepredicteq{1,iExp(n)},'y-.')
##%     plot(forcepredicteq{2,iExp(n)},'r-.')
##%     plot(forcepredicteq{3,iExp(n)},'c-.')
##%     plot(forcepredicteq{4,iExp(n)},'m--')
##%     plot(forcepredicteq{5,iExp(n)},'g--')
##%     plot(forcepredicteq{6,iExp(n)},'b--')
##
##    ylim(ylimit)
##    set(gca,'box','on','layer','top','Xtick',[])
##%     if any(n==[1,3,5]),ylabel('Force (N)'),end
##end
##width = 5.8;height = 2.4;
##% set(gcf,'units','inches','papersize',[width height],'position',[1 9 width height],'paperpositionmode','auto')
##% saveFolder = '/Users/charlesrichter/Desktop/CR Shared Drive/Thesis/figures/';
##% saveas(gcf,[saveFolder,'forces_usquared'],'pdf')
