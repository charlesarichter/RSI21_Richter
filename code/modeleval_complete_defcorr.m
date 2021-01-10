%% Model Evaluation
clear,clc,close all

pkg load optim
pkg load signal

##addpath("/home/crichter/ccsl/Code_CCSL_Drive/Code/Model Evaluation")
##addpath("/home/crichter/ccsl/NSR20_Richter/code/models")
addpath("/home/crichter/ccsl/RSI21_Richter/code/models")

%% Load Data
% folderName = '\\zeus\shared\Charlie Richter\FUTEK Experiments\Experiment 7-12\';
% [spreadsheet,~,~] = xlsread([folderName 'selecteddata_defcorr_all_70.xls']);
%----
folderName = '/home/crichter/ccsl/';
spreadsheet = dlmread('/home/crichter/ccsl/selecteddata_defcorr_all_70.csv',',',2,0);
%----
%% Experiments for Testing

exps = [8 26 33 37 39 48 52 53 95 103]; % 1,78 removed
% exps = [22,32,38,42,51,56,58,61,85,105]; % Make sure this is consistent with force_simulation.m

% totalexps = unique(spreadsheet(:,1));
% exps = totalexps(totalexps>0);

nExps = length(exps);

forcepredict = cell(1,nExps);
forcepredicteq = cell(1,nExps);
mae1 = nan(1,nExps);mae1eq = nan(1,nExps);
rmd1 = nan(1,nExps);rmd1eq = nan(1,nExps);
eqdata = cell(nExps,1);
force = cell(nExps,1);

% Initialize Coefficients - MEAN VALUES CALCULATED FROM *ALL* EXPS
beta_jw = nan(nExps,3); % [Cl Crotlift Cm22] = [1.6864 1.4752 0.7536], STD = [0.2389 0.9920 0.4233]
beta_rw_liftam = nan(nExps,3); % [Cl Camtrans Camrotddot]  = [1.8429 0.4382 -0.3152], STD = [0.2430 0.5158 1.0283]
beta_trans_rot = nan(nExps,2); % [Cl Crotlift] = [1.7361 -0.0908], STD = [0.2453 0.8649]
beta_trans_amrot = nan(nExps,2); % [Cl Camrot] = [1.7420 0.0755], STD = [0.2432 0.6126]
beta_trans_amtrans = nan(nExps,2); % [Cl Camtrans] = [1.7354 0.3756], STD = [0.2514 0.5499]
beta_trans_amtrans_amrot = nan(nExps,3); % [Cl Camtrans Camrot] = [1.7109 0.6684 -0.3172], STD = [0.2316 0.4530 0.5491]
beta_trans_amtransRW = nan(nExps,2);
beta_trans_amrotddot = nan(nExps,2); % [Cl Camrotddot] = [1.8360 -0.2577], STD = [0.2404 0.9987]

% EQ coefficients
beta_eq1nf = [1.5723863e-8 1.2146853e-7 1.2146853e-7 1.2146853e-7 7.0700918e-7];
beta_eq2nf = [1.419276e-8 1.2280232e-7 1.2280232e-7 1.1480129e-6];
beta_eq3nf = [1.4180717e-8 1.2409225e-7 1.0660193e-6];
beta_eq4nf = [1.4526972e-8 1.2489714e-7 1.0872654e-6];
beta_eq5nf = [1.2524723e-7 5.8678208e-11 1.1539232e-6];
beta_eq6nf = [1.3552796e-8 1.2505514e-7 1.1357268e-6];

%% Loop Through Test Experiments - With Coefficient Fitting
for i = 1:nExps
    disp(i)
    eqdata{i} = spreadsheet(spreadsheet(:,1)==exps(i),2:end);
    force{i} = eqdata{i}(:,1);

    % Planar Inertial
    forcepredict{1,i} = fn_planarinertial(eqdata{i});
    mae1(1,i) = mean(abs(force{i}-forcepredict{1,i}));rmd1(1,i) = mae1(1,i)/mean(force{i});

    % Total Inertial
    forcepredict{2,i} = fn_totalinertial(eqdata{i});
    mae1(2,i) = mean(abs(force{i}-forcepredict{2,i}));rmd1(2,i) = mae1(2,i)/mean(force{i});

    % Planar Lift Planar Inertial
    beta0_plpi = 1.5;
    beta_plpi(i) = nlinfit(eqdata{i},force{i},@fn_planarlift_planarinertial,beta0_plpi); 
    forcepredict{3,i} = fn_planarlift_planarinertial(beta_plpi(i),eqdata{i});
    mae1(3,i) = mean(abs(force{i}-forcepredict{3,i}));rmd1(3,i) = mae1(3,i)/mean(force{i});

    % Planar Lift Total Inertial
    beta0_plti = 1.5;
    beta_plti(i) = nlinfit(eqdata{i},force{i},@fn_planarlift_totalinertial,beta0_plti); 
    forcepredict{4,i} = fn_planarlift_totalinertial(beta_plti(i),eqdata{i});
    mae1(4,i) = mean(abs(force{i}-forcepredict{4,i}));rmd1(4,i) = mae1(4,i)/mean(force{i});

    % RW Lift Only Total Inertial
    beta0_rw_liftonly = 1.5;
    beta_rw_liftonly(i) = nlinfit(eqdata{i},force{i},@fn_rw_liftonly,beta0_rw_liftonly); %#ok<*SAGROW>
    forcepredict{5,i} = fn_rw_liftonly(beta_rw_liftonly(i),eqdata{i});
    mae1(5,i) = mean(abs(force{i}-forcepredict{5,i}));rmd1(5,i) = mae1(5,i)/mean(force{i});

    % RW Lift and AM Total Inertial
    beta0_rw_liftam = [1.5 1 1];
    beta_rw_liftam(i,:) = nlinfit(eqdata{i},force{i},@fn_rw_liftam,beta0_rw_liftam);
    forcepredict{6,i} = fn_rw_liftam(beta_rw_liftam(i,:),eqdata{i});
    mae1(6,i) = mean(abs(force{i}-forcepredict{6,i}));rmd1(6,i) = mae1(6,i)/mean(force{i});

    % JW Lift and Total Inertial
    beta0_jw = [1.5 pi 1];
    beta_jw(i,:) = nlinfit(eqdata{i},force{i},@fn_jw,beta0_jw);
    forcepredict{7,i} = fn_jw(beta_jw(i,:),eqdata{i});
    mae1(7,i) = mean(abs(force{i}-forcepredict{7,i}));rmd1(7,i) = mae1(7,i)/mean(force{i});

    % Analytical: fn_trans_rot -- COMPLEXITY 118
    beta0_trans_rot = [1.5 2];
    beta_trans_rot(i,:) = nlinfit(eqdata{i},force{i},@fn_trans_rot,beta0_trans_rot); %#ok<*SAGROW>
    forcepredict{8,i} = fn_trans_rot(beta_trans_rot(i,:),eqdata{i});
    mae1(8,i) = mean(abs(force{i}-forcepredict{8,i}));rmd1(8,i) = mae1(8,i)/mean(force{i});
    corrholder = corrcoef(force{i},forcepredict{8,i});corr1(8,i) = corrholder(2);

    % Analytical: fn_trans_amrot -- COMPLEXITY 113
    beta0_trans_amrot = [1.5 1];
    beta_trans_amrot(i,:) = nlinfit(eqdata{i},force{i},@fn_trans_amrot,beta0_trans_amrot); %#ok<*SAGROW>
    forcepredict{9,i} = fn_trans_amrot(beta_trans_amrot(i,:),eqdata{i});
    mae1(9,i) = mean(abs(force{i}-forcepredict{9,i}));rmd1(9,i) = mae1(9,i)/mean(force{i});
    corrholder = corrcoef(force{i},forcepredict{9,i});corr1(9,i) = corrholder(2);

    % Analytical: fn_trans_amtrans -- COMPLEXITY 115
    beta0_trans_amtrans = [1.5 1];
    beta_trans_amtrans(i,:) = nlinfit(eqdata{i},force{i},@fn_trans_amtrans,beta0_trans_amtrans); %#ok<*SAGROW>
    forcepredict{10,i} = fn_trans_amtrans(beta_trans_amtrans(i,:),eqdata{i});
    mae1(10,i) = mean(abs(force{i}-forcepredict{10,i}));rmd1(10,i) = mae1(10,i)/mean(force{i});
    corrholder = corrcoef(force{i},forcepredict{10,i});corr1(10,i) = corrholder(2);

    % Analytical: fn_trans_amtrans_amrot -- COMPLEXITY 120
    beta0_trans_amtrans_amrot = [1.5 1 1];
    beta_trans_amtrans_amrot(i,:) = nlinfit(eqdata{i},force{i},@fn_trans_amtrans_amrot,beta0_trans_amtrans_amrot); %#ok<*SAGROW>
    forcepredict{11,i} = fn_trans_amtrans_amrot(beta_trans_amtrans_amrot(i,:),eqdata{i});
    mae1(11,i) = mean(abs(force{i}-forcepredict{11,i}));rmd1(11,i) = mae1(11,i)/mean(force{i});
    corrholder = corrcoef(force{i},forcepredict{11,i});corr1(11,i) = corrholder(2);

    % Analytical: fn_trans_amrotddot -- COMPLEXITY 106
    beta0_trans_amrotddot = [1.5 1];
    beta_trans_amrotddot(i,:) = nlinfit(eqdata{i},force{i},@fn_trans_amrotddot,beta0_trans_amrotddot); %#ok<*SAGROW>
    forcepredict{12,i} = fn_trans_amrotddot(beta_trans_amrotddot(i,:),eqdata{i});
    mae1(12,i) = mean(abs(force{i}-forcepredict{12,i}));rmd1(12,i) = mae1(12,i)/mean(force{i});
    corrholder = corrcoef(force{i},forcepredict{12,i});corr1(12,i) = corrholder(2);

    % EQ 1
    beta0_eq1 = [1e-7 1e-7 1e-7 1e-7 1e-7];
    beta_eq1(i,:) = nlinfit(eqdata{i},force{i},@fn_eq1,beta0_eq1);
    forcepredicteq{1,i} = fn_eq1(beta_eq1(i,:),eqdata{i});
    mae1eq(1,i) = mean(abs(force{i}-forcepredicteq{1,i}));rmd1eq(1,i) = mae1eq(1,i)/mean(force{i});

    % EQ 2
    beta0_eq2 = [1e-7 1e-7 1e-7 1e-7];
    beta_eq2(i,:) = nlinfit(eqdata{i},force{i},@fn_eq2,beta0_eq2);
    forcepredicteq{2,i} = fn_eq2(beta_eq2(i,:),eqdata{i});
    mae1eq(2,i) = mean(abs(force{i}-forcepredicteq{2,i}));rmd1eq(2,i) = mae1eq(2,i)/mean(force{i});

    % EQ 3
    beta0_eq3 = [1e-7 1e-7 1e-7];
    beta_eq3(i,:) = nlinfit(eqdata{i},force{i},@fn_eq3,beta0_eq3);
    forcepredicteq{3,i} = fn_eq3(beta_eq3(i,:),eqdata{i});
    mae1eq(3,i) = mean(abs(force{i}-forcepredicteq{3,i}));rmd1eq(3,i) = mae1eq(3,i)/mean(force{i});
    
    % EQ 4
    beta0_eq4 = [1e-7 1e-7 1e-7];
    beta_eq4(i,:) = nlinfit(eqdata{i},force{i},@fn_eq4,beta0_eq4);
    forcepredicteq{4,i} = fn_eq4(beta_eq4(i,:),eqdata{i});
    mae1eq(4,i) = mean(abs(force{i}-forcepredicteq{4,i}));rmd1eq(4,i) = mae1eq(4,i)/mean(force{i});
    
    % EQ 5
    beta0_eq5 = [1e-7 1e-7 1e-7];
    beta_eq5(i,:) = nlinfit(eqdata{i},force{i},@fn_eq5,beta0_eq5);
    forcepredicteq{5,i} = fn_eq5(beta_eq5(i,:),eqdata{i});
    mae1eq(5,i) = mean(abs(force{i}-forcepredicteq{5,i}));rmd1eq(5,i) = mae1eq(5,i)/mean(force{i});
    
    % EQ 6
    beta0_eq6 = [1e-7 1e-7 1e-7];
    beta_eq6(i,:) = nlinfit(eqdata{i},force{i},@fn_eq6,beta0_eq6);
    forcepredicteq{6,i} = fn_eq6(beta_eq6(i,:),eqdata{i});
    mae1eq(6,i) = mean(abs(force{i}-forcepredicteq{6,i}));rmd1eq(6,i) = mae1eq(6,i)/mean(force{i});
end
%% Averages
fitmae_an = mean(mae1,2);
fitrmd_an = mean(rmd1,2);
fitmae_eq = mean(mae1eq,2);
fitrmd_eq = mean(rmd1eq,2);

%% Coeff. Averages
beta_plpi_nf = mean(beta_plpi);
beta_plti_nf = mean(beta_plti);
beta_rw_liftonly_nf = mean(beta_rw_liftonly);

beta_rw_liftam_nf = mean(beta_rw_liftam,1);
beta_rw_std = std(beta_rw_liftam,1);
beta_jw_nf = mean(beta_jw,1);
beta_jw_std = std(beta_jw,1);

beta_trans_rot_nf = mean(beta_trans_rot,1);
beta_trans_amrot_nf = mean(beta_trans_amrot,1);
beta_trans_amtrans_nf = mean(beta_trans_amtrans,1);
beta_trans_amtrans_amrot_nf = mean(beta_trans_amtrans_amrot,1);
beta_trans_amrotddot_nf = mean(beta_trans_amrotddot,1);

clear mae1 mae1eq rmd1 rmd1eq eqdata force

%% Loop Through Test Experiments - No Coefficient Fitting
forcepredict2 = cell(1,nExps);
forcepredicteq2 = cell(1,nExps);
mae2 = nan(1,nExps);mae2eq = nan(1,nExps);
rmd2 = nan(1,nExps);rmd2eq = nan(1,nExps);
eqdata = cell(nExps,1);
force = cell(nExps,1);

for i = 1:nExps
    eqdata{i} = spreadsheet(spreadsheet(:,1)==exps(i),2:end);
    force{i} = eqdata{i}(:,1);

    % Planar Inertial -- COMPLEXITY 13
    forcepredict2{1,i} = fn_planarinertial(eqdata{i});
    mae2(1,i) = mean(abs(force{i}-forcepredict2{1,i}));rmd2(1,i) = mae2(1,i)/mean(force{i});
    
    % Total Inertial -- COMPLEXITY 27
    forcepredict2{2,i} = fn_totalinertial(eqdata{i});
    mae2(2,i) = mean(abs(force{i}-forcepredict2{2,i}));rmd2(2,i) = mae2(2,i)/mean(force{i});

    % Planar Lift Planar Inertial -- COMPLEXITY 32
%     beta_plpi_nf = 1.7942;
    forcepredict2{3,i} = fn_planarlift_planarinertial(beta_plpi_nf,eqdata{i});
    mae2(3,i) = mean(abs(force{i}-forcepredict2{3,i}));rmd2(3,i) = mae2(3,i)/mean(force{i});

    % Planar Lift Total Inertial -- COMPLEXITY 46
%     beta_plti_nf = 1.7327;
    forcepredict2{4,i} = fn_planarlift_totalinertial(beta_plti_nf,eqdata{i});
    mae2(4,i) = mean(abs(force{i}-forcepredict2{4,i}));rmd2(4,i) = mae2(4,i)/mean(force{i});

    % RW Lift Only Total Inertial -- COMPLEXITY 96
%     beta_rw_liftonly_nf = 1.7410;
    forcepredict2{5,i} = fn_rw_liftonly(beta_rw_liftonly_nf,eqdata{i});
    mae2(5,i) = mean(abs(force{i}-forcepredict2{5,i}));rmd2(5,i) = mae2(5,i)/mean(force{i});

    % RW Lift and AM Total Inertial -- COMPLEXITY 121
%     beta_rw_liftam_nf = [1.8429 0.4382 -0.3152];
    forcepredict2{6,i} = fn_rw_liftam(beta_rw_liftam_nf,eqdata{i});
    mae2(6,i) = mean(abs(force{i}-forcepredict2{6,i}));rmd2(6,i) = mae2(6,i)/mean(force{i});

    % JW Lift and Total Inertial -- COMPLEXITY 130
%     beta_jw_nf = [1.6864 1.4752 0.7536];
    forcepredict2{7,i} = fn_jw(beta_jw_nf,eqdata{i});
    mae2(7,i) = mean(abs(force{i}-forcepredict2{7,i}));rmd2(7,i) = mae2(7,i)/mean(force{i});

    % Analytical: fn_trans_rot -- COMPLEXITY 118
%     beta_trans_rot_nf = [1.7361 -0.0908];
    forcepredict2{8,i} = fn_trans_rot(beta_trans_rot_nf,eqdata{i});
    mae2(8,i) = mean(abs(force{i}-forcepredict2{8,i}));rmd2(8,i) = mae2(8,i)/mean(force{i});

    % Analytical: fn_trans_amrot -- COMPLEXITY 113
%     beta_trans_amrot_nf = [1.7361 -0.0908];
    forcepredict2{9,i} = fn_trans_amrot(beta_trans_amrot_nf,eqdata{i});
    mae2(9,i) = mean(abs(force{i}-forcepredict2{9,i}));rmd2(9,i) = mae2(9,i)/mean(force{i});

    % Analytical: fn_trans_amtrans -- COMPLEXITY 115
%     beta_trans_amtrans_nf = [1.7420 0.0755];
    forcepredict2{10,i} = fn_trans_amtrans(beta_trans_amtrans_nf,eqdata{i});
    mae2(10,i) = mean(abs(force{i}-forcepredict2{10,i}));rmd2(10,i) = mae2(10,i)/mean(force{i});

    % Analytical: fn_trans_amtrans_amrot -- COMPLEXITY 120
%     beta_trans_amtrans_amrot_nf = [1.7109 0.6684 -0.3172];
    forcepredict2{11,i} = fn_trans_amtrans_amrot(beta_trans_amtrans_amrot_nf,eqdata{i});
    mae2(11,i) = mean(abs(force{i}-forcepredict2{11,i}));rmd2(11,i) = mae2(11,i)/mean(force{i});

    % Analytical: fn_trans_amrotddot -- COMPLEXITY 106
%     beta_trans_amrotddot_nf = [1.8360 -0.2577];
    forcepredict2{12,i} = fn_trans_amrotddot(beta_trans_amrotddot_nf,eqdata{i});
    mae2(12,i) = mean(abs(force{i}-forcepredict2{12,i}));rmd2(12,i) = mae2(12,i)/mean(force{i});

    % EQ 1 -- COMPLEXITY 24
    forcepredicteq2{1,i} = fn_eq1(beta_eq1nf,eqdata{i});
    mae2eq(1,i) = mean(abs(force{i}-forcepredicteq2{1,i}));rmd2eq(1,i) = mae2eq(1,i)/mean(force{i});
    
    % EQ 2 -- COMPLEXITY 19
    forcepredicteq2{2,i} = fn_eq2(beta_eq2nf,eqdata{i});
    mae2eq(2,i) = mean(abs(force{i}-forcepredicteq2{2,i}));rmd2eq(2,i) = mae2eq(2,i)/mean(force{i});
    
    % EQ 3 -- COMPLEXITY 16
    forcepredicteq2{3,i} = fn_eq3(beta_eq3nf,eqdata{i});
    mae2eq(3,i) = mean(abs(force{i}-forcepredicteq2{3,i}));rmd2eq(3,i) = mae2eq(3,i)/mean(force{i});
    
    % EQ 4 -- COMPLEXITY 14
    forcepredicteq2{4,i} = fn_eq4(beta_eq4nf,eqdata{i});
    mae2eq(4,i) = mean(abs(force{i}-forcepredicteq2{4,i}));rmd2eq(4,i) = mae2eq(4,i)/mean(force{i});
    
    % EQ 5 -- COMPLEXITY 13
    forcepredicteq2{5,i} = fn_eq5(beta_eq5nf,eqdata{i});
    mae2eq(5,i) = mean(abs(force{i}-forcepredicteq2{5,i}));rmd2eq(5,i) = mae2eq(5,i)/mean(force{i});
    
    % EQ 6 -- COMPLEXITY 12
    forcepredicteq2{6,i} = fn_eq6(beta_eq6nf,eqdata{i});
    mae2eq(6,i) = mean(abs(force{i}-forcepredicteq2{6,i}));rmd2eq(6,i) = mae2eq(6,i)/mean(force{i});
end
%% Averages
nofitmae_an = mean(mae2,2);
nofitrmd_an = mean(rmd2,2);
nofitmae_eq = mean(mae2eq,2);
nofitrmd_eq = mean(rmd2eq,2);

%% Polynomial
pOrder = 2:1:30;
for n = 1:length(pOrder)
    for k = 1:nExps
        [pks,locs]=findpeaks(abs(eqdata{k}(:,2)),'minpeakheight',0);
        middle = round(length(eqdata{k}(:,2))/2);        
        [middleval,middlepos] = min(abs(locs-middle));
        fpoly{k} = force{k}(locs(middlepos-1):locs(middlepos+1));
        
        polyx{k} = (1:length(fpoly{k}))';
        [p{k},s{k},mu{k}] = polyfit(polyx{k},fpoly{k},pOrder(n));
        y{k} = polyval(p{k},polyx{k},[],mu{k});
        e{k} = fpoly{k}-y{k};
        maepoly{k} = mean(abs(e{k}));
    end
    polymae(n) = mean(cell2mat(maepoly));

%     figure(4)
%     for m = 1:nExps
%         subplot(nExps,1,m),hold on
%         plot(polyx{m},fpoly{m},'k','linewidth',2)
%         plot(polyx{m},y{m},'b')
% %         plot(polyx{m},e{m},'g')
%         ylim([-.015 .015]),hold off
%     end
%     figure(3),hold on
%     plot(3*n-1,polymae(n),'k+')
end
% plot(3*pOrder-1,polymae,'k','linewidth',2)

%% Complexity
complex_an = [13 27 32 46 96 121 130 118 113 115 120 106];
complex_eq = [24 19 16 14 13 12];

%% Pareto - Fitted Coefficients
figure(1),clf,hold on
% set(gcf,'units','inches','outerposition',[0 0 5.8 3])
smallmsize = 16;
rotangle = 30;
fsize = 6;

%==== FOR PNAS - SWITCHED FROM NOFITMAE TO FITMAE
% plot(complex_an(3),fitmae_an(3),'b.','Markersize',20)%smallmsize)
% text(complex_an(3),fitmae_an(3),'   U^2 Planar 1','rotation',rotangle,'fontsize',fsize)
% plot(complex_an(4),fitmae_an(4),'g.','Markersize',20)%,smallmsize)
% text(complex_an(4),fitmae_an(4),'   U^2 Planar 2','rotation',rotangle,'fontsize',fsize)
% 
% plot(complex_an(5),fitmae_an(5),'c.','Markersize',20)%,smallmsize)
% text(complex_an(5),fitmae_an(5),'   U^2 Total','rotation',rotangle,'fontsize',fsize)
% 
% % plot(complex_an(8:12),fitmae_an(8:12),'g.','Markersize',smallmsize)
% % text(complex_an(11)-4,fitmae_an(11)+4e-5,'   Other Combinations','rotation',rotangle,'fontsize',fsize)
% 
% plot(complex_an(6),fitmae_an(6),'r.','Markersize',20)
% text(complex_an(6),fitmae_an(6)+2e-5,'   Whitney-Wood','rotation',rotangle,'fontsize',fsize)
% 
% plot(complex_an(7),fitmae_an(7),'m.','Markersize',20)
% text(complex_an(7),fitmae_an(7),'   Pesavento-Wang','rotation',rotangle,'fontsize',fsize)
% 
% plot(complex_eq([1,2,3,4,6]),fitmae_eq([1,2,3,4,6]),'k.','Markersize',smallmsize)
% text(complex_eq(1)-6,fitmae_eq(1)-1.5e-4,{'EQ';'Models'},'horizontalalignment','center','fontsize',fsize)
%====

%==== 
plot(complex_an(3),nofitmae_an(3),'b.','Markersize',20)%smallmsize)
text(complex_an(3),nofitmae_an(3),'   U^2 Planar 1','rotation',rotangle,'fontsize',fsize)
plot(complex_an(4),nofitmae_an(4),'g.','Markersize',20)%,smallmsize)
text(complex_an(4),nofitmae_an(4),'   U^2 Planar 2','rotation',rotangle,'fontsize',fsize)

plot(complex_an(5),nofitmae_an(5),'c.','Markersize',20)%,smallmsize)
text(complex_an(5),nofitmae_an(5),'   U^2 Total','rotation',rotangle,'fontsize',fsize)

% plot(complex_an(8:12),nofitmae_an(8:12),'g.','Markersize',smallmsize)
% text(complex_an(11)-4,nofitmae_an(11)+4e-5,'   Other Combinations','rotation',rotangle,'fontsize',fsize)

plot(complex_an(6),nofitmae_an(6),'r.','Markersize',20)
text(complex_an(6),nofitmae_an(6)+2e-5,'   Whitney-Wood','rotation',rotangle,'fontsize',fsize)

plot(complex_an(7),nofitmae_an(7),'m.','Markersize',20)
text(complex_an(7),nofitmae_an(7),'   Pesavento-Wang','rotation',rotangle,'fontsize',fsize)

plot(complex_eq([1,2,3,4,6]),nofitmae_eq([1,2,3,4,6]),'k.','Markersize',smallmsize)
text(complex_eq(1)-6,nofitmae_eq(1)-1.5e-4,{'EQ';'Models'},'horizontalalignment','center','fontsize',fsize)

% plot(complex_an(3:end),fitmae_an(3:end),'ro')
% plot(complex_eq,fitmae_eq,'bo')
% % % plot(3*pOrder-1,polymae,'k.','markersize',6)
%====

ylim([0 19e-4])
xlim([0 180])
ylabel('Mean Absolute Error (10^-^3 N)','fontsize',fsize)
xlabel('Equation Size (number of operators)','fontsize',fsize)
width=8.7;height=5;
lbuf=1;rbuf=.2;
bbuf=.8;tbuf=0;
set(gca,'box','on','layer','top','fontunits','points','fontsize',fsize)
set(gca,'units','centimeters','position',[lbuf bbuf width-lbuf-rbuf height-bbuf-tbuf])
set(gcf,'units','centimeters','paperunits','centimeters','papersize',[width height],...
    'paperposition',[0 0 width height])

% saveFolder = '/Users/charlesrichter/Desktop/CR Shared Drive/Thesis/figures/';
##saveFolder = '/Users/charlesrichter/Desktop/PNAS11_Richter/figures/';
##saveFolder = '/home/crichter/ccsl/NSR20_Richter/paper/figures/';
saveFolder = '/home/crichter/ccsl/RSI21_Richter/paper/figures/';
saveas(gcf,[saveFolder,'pareto_mae_nofit_PNAS'],'epsc')
% MARKERSIZES, FONTSIZES, DIMENSIONS, ETC ARE CHANGED FROM THE ORIGINAL
% THESIS PLOT TO BE USED IN THE JOURNAL PUBLICATION.  REFER TO THE PLOTS
% BELOW TO RETURN TO THE ORIGINAL VALUES.


%% Pareto - Fitted Coefficients
figure(2),clf,hold on
% set(gcf,'units','inches','outerposition',[0 0 5.8 3])
smallmsize = 16;
rotangle = 30;
fsize = 8;

plot(complex_an(3),fitmae_an(3),'b.','Markersize',smallmsize)
text(complex_an(3),fitmae_an(3),'  U^2 Planar 1','rotation',rotangle,'fontsize',fsize)
plot(complex_an(4),fitmae_an(4),'b.','Markersize',smallmsize)
text(complex_an(4),fitmae_an(4),'  U^2 Planar 2','rotation',rotangle,'fontsize',fsize)

plot(complex_an(5),fitmae_an(5),'c.','Markersize',smallmsize)
text(complex_an(5),fitmae_an(5),'  U^2 Total','rotation',rotangle,'fontsize',fsize)

plot(complex_an(8:12),fitmae_an(8:12),'g.','Markersize',smallmsize)
text(complex_an(10),fitmae_an(10)+2e-5,' Other Combinations','rotation',rotangle,'fontsize',fsize)

plot(complex_an(6),fitmae_an(6),'r.','Markersize',24)
text(complex_an(6),fitmae_an(6),'   Whitney-Wood','rotation',rotangle,'fontsize',fsize)

plot(complex_an(7),fitmae_an(7),'m.','Markersize',24)
text(complex_an(7),fitmae_an(7),'   Pesavento-Wang','rotation',rotangle,'fontsize',fsize)

plot(complex_eq,fitmae_eq,'k.','Markersize',smallmsize)
text(complex_eq(1)-6,fitmae_eq(1)-1.2e-4,{'EQ';'Models'},'horizontalalignment','center','fontsize',fsize)

% plot(complex_an(3:end),fitmae_an(3:end),'ro')
% plot(complex_eq,fitmae_eq,'bo')
plot(3*pOrder-1,polymae,'k.','markersize',6)

ylim([0 17e-4])
xlim([0 170])
ylabel('Mean Absolute Error')
xlabel('Equation Size (number of operators)')
width=5.8;height=3;
    set(gca,'box','on','layer','top')
set(gcf,'units','inches','papersize',[width height],'position',[1 9 width height],'paperpositionmode','auto')

% saveFolder = '/Users/charlesrichter/Desktop/CR Shared Drive/Thesis/figures/';
% saveas(gcf,[saveFolder,'pareto_mae_oct23'],'pdf')


%% Force Plots
figure(3),clf
set(gcf,'defaultlinelinewidth',.5)
iExp = [2 4 6 9];
shadecolor = .9*[1,1,1];

for n = 1:length(iExp)
    subplot(2,2,n),hold on  
    strokedot = eqdata{iExp(n)}(:,3);
    cross1 = find([0;diff(sign(strokedot))]);
    if cross1(1)<100
        fill([cross1(1) cross1(2) cross1(2) cross1(1)],100*[-1 -1 1 1],shadecolor,'edgecolor',shadecolor)
        fill([cross1(3) cross1(4) cross1(4) cross1(3)],100*[-1 -1 1 1],shadecolor,'edgecolor',shadecolor)
        xlim([cross1(1) length(force{iExp(n)})])
    else
        fill([0 cross1(1) cross1(1) 0],100*[-1 -1 1 1],shadecolor,'edgecolor',shadecolor)
        fill([cross1(2) cross1(3) cross1(3) cross1(2)],100*[-1 -1 1 1],shadecolor,'edgecolor',shadecolor)
        xlim([1 length(force{iExp(n)})])
    end
    clear strokedot
    plot(force{iExp(n)},'k')
    ylimit = 1.2*[min(force{iExp(n)}) max(force{iExp(n)})];
%     ylimit = 1.2*[min(force{iExp(n)}) max(forcepredict{6,iExp(n)})];

    plot(forcepredict{3,iExp(n)},'b')
    plot(forcepredict{4,iExp(n)},'r')

%     plot(forcepredict{5,iExp(n)},'c')
%     plot(forcepredict{6,iExp(n)},'m')

%     plot(forcepredict{7,iExp(n)},'r')

%     plot(forcepredicteq{1,iExp(n)},'y-.')
%     plot(forcepredicteq{2,iExp(n)},'r-.')
%     plot(forcepredicteq{3,iExp(n)},'c-.')
%     plot(forcepredicteq{4,iExp(n)},'m--')
%     plot(forcepredicteq{5,iExp(n)},'g--')
%     plot(forcepredicteq{6,iExp(n)},'b--')

    ylim(ylimit)
    set(gca,'box','on','layer','top','Xtick',[])
%     if any(n==[1,3,5]),ylabel('Force (N)'),end
end
width = 5.8;height = 2.4;
% set(gcf,'units','inches','papersize',[width height],'position',[1 9 width height],'paperpositionmode','auto')
% saveFolder = '/Users/charlesrichter/Desktop/CR Shared Drive/Thesis/figures/';
% saveas(gcf,[saveFolder,'forces_usquared'],'pdf')