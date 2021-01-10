clear,clc,close all

pkg load optim
pkg load signal

addpath("/home/crichter/ccsl/RSI21_Richter/code/models")

%% Load Data
folderName = '/home/crichter/ccsl/';
spreadsheet = dlmread('/home/crichter/ccsl/selecteddata_defcorr_all_70.csv',',',2,0);

%% Determine training, testing and other removed experiment IDs
exps_all = unique(spreadsheet(:,1));
exps_testing = [8 26 33 37 39 48 52 53 95 103]; % 1,78 removed
exps_removed = [1, 78];
exps_training = setdiff(exps_all, [exps_testing, exps_removed, 0]);

% Mask off data to only include rows for which the experiment number (first column) is not in exps_testing, exps_removed, or 0.
should_include = any(repmat(spreadsheet(:,1),1,size(exps_training)(1)) == repmat(exps_training',size(spreadsheet)(1),1),2);
size(should_include)

eqdata = spreadsheet(should_include,2:end);
force = spreadsheet(should_include,2);

##forcepredict = cell(1,nExps);
##forcepredicteq = cell(1,nExps);
##mae1 = nan(1,nExps);mae1eq = nan(1,nExps);
##rmd1 = nan(1,nExps);rmd1eq = nan(1,nExps);
##eqdata = cell(nExps,1);
##force = cell(nExps,1);

% Initialize Coefficients - MEAN VALUES CALCULATED FROM *ALL* EXPS
##beta_jw = nan(nExps,3); % [Cl Crotlift Cm22] = [1.6864 1.4752 0.7536], STD = [0.2389 0.9920 0.4233]
##beta_rw_liftam = nan(nExps,3); % [Cl Camtrans Camrotddot]  = [1.8429 0.4382 -0.3152], STD = [0.2430 0.5158 1.0283]
##beta_trans_rot = nan(nExps,2); % [Cl Crotlift] = [1.7361 -0.0908], STD = [0.2453 0.8649]
##beta_trans_amrot = nan(nExps,2); % [Cl Camrot] = [1.7420 0.0755], STD = [0.2432 0.6126]
##beta_trans_amtrans = nan(nExps,2); % [Cl Camtrans] = [1.7354 0.3756], STD = [0.2514 0.5499]
##beta_trans_amtrans_amrot = nan(nExps,3); % [Cl Camtrans Camrot] = [1.7109 0.6684 -0.3172], STD = [0.2316 0.4530 0.5491]
##beta_trans_amtransRW = nan(nExps,2);
##beta_trans_amrotddot = nan(nExps,2); % [Cl Camrotddot] = [1.8360 -0.2577], STD = [0.2404 0.9987]

%% Attempt coefficient fitting using all training datasets simultaneously

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
beta0_plpi = 1.5;
beta_plpi = nlinfit(eqdata,force,@fn_planarlift_planarinertial,beta0_plpi)
##    forcepredict{3,i} = fn_planarlift_planarinertial(beta_plpi(i),eqdata{i});
##    mae1(3,i) = mean(abs(force{i}-forcepredict{3,i}));rmd1(3,i) = mae1(3,i)/mean(force{i});
##
##    % Planar Lift Total Inertial
##    beta0_plti = 1.5;
##    beta_plti(i) = nlinfit(eqdata{i},force{i},@fn_planarlift_totalinertial,beta0_plti); 
beta0_plti = 1.5;
beta_plti = nlinfit(eqdata,force,@fn_planarlift_totalinertial,beta0_plti)
##    forcepredict{4,i} = fn_planarlift_totalinertial(beta_plti(i),eqdata{i});
##    mae1(4,i) = mean(abs(force{i}-forcepredict{4,i}));rmd1(4,i) = mae1(4,i)/mean(force{i});
##
##    % RW Lift Only Total Inertial
##    beta0_rw_liftonly = 1.5;
##    beta_rw_liftonly(i) = nlinfit(eqdata{i},force{i},@fn_rw_liftonly,beta0_rw_liftonly); %#ok<*SAGROW>
beta0_rw_liftonly = 1.5;
beta_rw_liftonly = nlinfit(eqdata,force,@fn_rw_liftonly,beta0_rw_liftonly)
##    forcepredict{5,i} = fn_rw_liftonly(beta_rw_liftonly(i),eqdata{i});
##    mae1(5,i) = mean(abs(force{i}-forcepredict{5,i}));rmd1(5,i) = mae1(5,i)/mean(force{i});
##
##    % RW Lift and AM Total Inertial
##    beta0_rw_liftam = [1.5 1 1];
##    beta_rw_liftam(i,:) = nlinfit(eqdata{i},force{i},@fn_rw_liftam,beta0_rw_liftam);
beta0_rw_liftam = [1.5 1 1];
beta_rw_liftam = nlinfit(eqdata,force,@fn_rw_liftam,beta0_rw_liftam)
##    forcepredict{6,i} = fn_rw_liftam(beta_rw_liftam(i,:),eqdata{i});
##    mae1(6,i) = mean(abs(force{i}-forcepredict{6,i}));rmd1(6,i) = mae1(6,i)/mean(force{i});
##
##    % JW Lift and Total Inertial
##    beta0_jw = [1.5 pi 1];
##    beta_jw(i,:) = nlinfit(eqdata{i},force{i},@fn_jw,beta0_jw);
beta0_jw = [1.5 pi 1];
beta_jw = nlinfit(eqdata,force,@fn_jw,beta0_jw)
##    forcepredict{7,i} = fn_jw(beta_jw(i,:),eqdata{i});
##    mae1(7,i) = mean(abs(force{i}-forcepredict{7,i}));rmd1(7,i) = mae1(7,i)/mean(force{i});
##
##    % Analytical: fn_trans_rot -- COMPLEXITY 118
##    beta0_trans_rot = [1.5 2];
##    beta_trans_rot(i,:) = nlinfit(eqdata{i},force{i},@fn_trans_rot,beta0_trans_rot); %#ok<*SAGROW>
beta0_trans_rot = [1.5 2];
beta_trans_rot = nlinfit(eqdata,force,@fn_trans_rot,beta0_trans_rot)
##    forcepredict{8,i} = fn_trans_rot(beta_trans_rot(i,:),eqdata{i});
##    mae1(8,i) = mean(abs(force{i}-forcepredict{8,i}));rmd1(8,i) = mae1(8,i)/mean(force{i});
##    corrholder = corrcoef(force{i},forcepredict{8,i});corr1(8,i) = corrholder(2);
##
##    % Analytical: fn_trans_amrot -- COMPLEXITY 113
##    beta0_trans_amrot = [1.5 1];
##    beta_trans_amrot(i,:) = nlinfit(eqdata{i},force{i},@fn_trans_amrot,beta0_trans_amrot); %#ok<*SAGROW>
beta0_trans_amrot = [1.5 1];
beta_trans_amrot = nlinfit(eqdata,force,@fn_trans_amrot,beta0_trans_amrot)
##    forcepredict{9,i} = fn_trans_amrot(beta_trans_amrot(i,:),eqdata{i});
##    mae1(9,i) = mean(abs(force{i}-forcepredict{9,i}));rmd1(9,i) = mae1(9,i)/mean(force{i});
##    corrholder = corrcoef(force{i},forcepredict{9,i});corr1(9,i) = corrholder(2);
##
##    % Analytical: fn_trans_amtrans -- COMPLEXITY 115
##    beta0_trans_amtrans = [1.5 1];
##    beta_trans_amtrans(i,:) = nlinfit(eqdata{i},force{i},@fn_trans_amtrans,beta0_trans_amtrans); %#ok<*SAGROW>
beta0_trans_amtrans = [1.5 1];
beta_trans_amtrans = nlinfit(eqdata,force,@fn_trans_amtrans,beta0_trans_amtrans)
##    forcepredict{10,i} = fn_trans_amtrans(beta_trans_amtrans(i,:),eqdata{i});
##    mae1(10,i) = mean(abs(force{i}-forcepredict{10,i}));rmd1(10,i) = mae1(10,i)/mean(force{i});
##    corrholder = corrcoef(force{i},forcepredict{10,i});corr1(10,i) = corrholder(2);
##
##    % Analytical: fn_trans_amtrans_amrot -- COMPLEXITY 120
##    beta0_trans_amtrans_amrot = [1.5 1 1];
##    beta_trans_amtrans_amrot(i,:) = nlinfit(eqdata{i},force{i},@fn_trans_amtrans_amrot,beta0_trans_amtrans_amrot); %#ok<*SAGROW>
beta0_trans_amtrans_amrot = [1.5 1 1];
beta_trans_amtrans_amrot = nlinfit(eqdata,force,@fn_trans_amtrans_amrot,beta0_trans_amtrans_amrot)
##    forcepredict{11,i} = fn_trans_amtrans_amrot(beta_trans_amtrans_amrot(i,:),eqdata{i});
##    mae1(11,i) = mean(abs(force{i}-forcepredict{11,i}));rmd1(11,i) = mae1(11,i)/mean(force{i});
##    corrholder = corrcoef(force{i},forcepredict{11,i});corr1(11,i) = corrholder(2);
##
##    % Analytical: fn_trans_amrotddot -- COMPLEXITY 106
##    beta0_trans_amrotddot = [1.5 1];
##    beta_trans_amrotddot(i,:) = nlinfit(eqdata{i},force{i},@fn_trans_amrotddot,beta0_trans_amrotddot); %#ok<*SAGROW>
beta0_trans_amrotddot = [1.5 1];
beta_trans_amrotddot = nlinfit(eqdata,force,@fn_trans_amrotddot,beta0_trans_amrotddot)
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
##clear mae1 mae1eq rmd1 rmd1eq eqdata force