%% Data Registration and Packaging Code for 3D Printed Wing Experiments
% Charles Richter
% Revised: July 10, 2011

% function [] = fnregister(folderName,toggleRegisterPlots,toggleEndPlots,kBandwidth)
clear,clc,close all,toggleRegisterPlots=1;toggleEndPlots=1;kBandwidth=7.5;
##folderName = '\\zeus\shared\Charlie Richter\FUTEK Experiments\Experiment 7-12\Wing 100sp 50lc 10uc for Thesis Plots\Notch 1\';
% folderName ='\\zeus\shared\Charlie Richter\FUTEK Experiments\Experiment 7-12\Wing 40sp 60lc 10uc\Notch 1\';
folderName = '/home/crichter/ccsl/Wing 100sp 50lc 10uc/Wing 100sp 50lc 10uc/Notch 1/';

pkg load signal

%% Import data
load ([folderName 'compdata.mat']);kinematics(isnan(kinematics)) = 0; %#ok<NODEF>
load ([folderName 'forcewindows.mat']);fwin = forcewindows;
load ([folderName 'massprops.mat'])
load ([folderName 'constants.mat'])
spanTotal = span+innerSpan;chordTotal = chordLower+chordUpper;
fwin = fwin-(impulse1force-1); % shift the ranges to post-truncation indices

%% Representative Force Plots for Thesis Showing Filtering and PSD
% flapfreq = 3.0518; % for force window #5
% window = fwin(5,1):fwin(5,2);
% data = force;
% 
% % Filter force
% cutoff = flapfreq*kBandwidth;width = 4;
% h=fdesign.lowpass('Fp,Fst,Ap,Ast',cutoff,cutoff+width,.1,60,1000);
% d=design(h,'equiripple'); %Lowpass FIR filter
% dataf1 = filtfilt(d.Numerator,1,data); %zero-phase filtering
% dataf = dataf1(window);
% 
% h2=fdesign.lowpass('Fp,Fst,Ap,Ast',50,54,.1,60,1000);
% d2=design(h2,'equiripple'); %Lowpass FIR filter
% dataf12 = filtfilt(d2.Numerator,1,data); %zero-phase filtering
% dataf2 = dataf12(window);
% 
% % Compute Amplitude Spectrum
% Fs = 1000;                    % Sampling frequency
% T = 1/Fs;                     % Sample time
% L = length(dataf);            % Length of signal
% t = (0:L-1)*T;                % Time vector
% y = data(window);
% 
% NFFT = 2^nextpow2(L); % Next power of 2 from length of y
% Y = fft(y,NFFT)/L;
% f = Fs/2*linspace(0,1,NFFT/2+1);
% 
% NFFT2 = 2^nextpow2(L); % Next power of 2 from length of y
% Y2 = fft(dataf,NFFT2)/L;
% f2 = Fs/2*linspace(0,1,NFFT2/2+1);
% 
% NFFT3 = 2^nextpow2(L); % Next power of 2 from length of y
% Y3 = fft(dataf2,NFFT3)/L;
% f3 = Fs/2*linspace(0,1,NFFT3/2+1);
% 
% %% Plot
% figure(1),clf,hold on
% plot(y,'c','linewidth',1)
% plot(dataf2,'b','linewidth',2)
% plot(dataf,'r','linewidth',4)
% plot(kinematics(window,1)*.06)
% % xlabel('Time (ms)')
% ylabel('Force (N)')
% hold off
% xlim([1 999])
% ylim([-.08 .08])
% % xlabel('Time (ms)')
% 
% ylabel('Measured Lift Force (N)')
% width = 5.8;height = 3;
% set(gca,'units','inches','position',[1 1 width height],'xtick',[])
% set(gcf,'units','inches','position',[-14 1 width+2 height+1.5])
% set(gcf,'papersize',[width+1.5 height+1],'paperpositionmode','auto')
% saveFolder = '\\zeus\shared\Charlie Richter\Thesis\figures\';
% % saveas(gcf,[saveFolder,'force1'],'png')
% 
% figure(2),clf,hold on
% plot(f,2*abs(Y(1:NFFT/2+1)),'c','linewidth',2),hold on
% plot(f3,2*abs(Y3(1:NFFT3/2+1)),'b','linewidth',1.5)
% plot(f2,2*abs(Y2(1:NFFT2/2+1)),'r','linewidth',1),hold off
% xlim([1 79])
% % ylim([0 .02])
% % title('Single-Sided Amplitude Spectrum of y(t)')
% xlabel('Frequency (Hz)')
% ylabel('Normalized Amplitude of FFT')
% set(gca,'units','inches','position',[1 1 width height],'XAxisLocation','top','box','on')
% set(gcf,'units','inches','position',[-14 1 width+2 height+1.5])
% set(gcf,'papersize',[width+1.5 height+1],'paperpositionmode','auto')
% % saveas(gcf,[saveFolder,'ampspec2'],'png')
% 

%% Obtain flapping frequency in domain of interest
% Fs = 1000;L = length(window); % sampling freq. and signal length
% NFFT = 2^nextpow2(L); % Next power of 2 from length of y
% f = Fs/2*linspace(0,1,NFFT/2+1); % Associated frequency vector
% 
%  % Lowpass filter force and kinematic data
% cutoff = flapfreq*kBandwidth;width = 4;
% h=fdesign.lowpass('Fp,Fst,Ap,Ast',cutoff,cutoff+width,.01,60,1000);
% d=design(h,'equiripple'); %Lowpass FIR filter
% forcef = filtfilt(d.Numerator,1,force); %zero-phase filtering
% 
% % Compute PSD of flapping in region
% Yforce = fft(force(window),NFFT)/L;psdforce = 2*abs(Yforce(1:NFFT/2+1));
% Yforcef = fft(forcef(window),NFFT)/L;psdforcef = 2*abs(Yforcef(1:NFFT/2+1));
% %%
% figure(1),clf,hold on
% plot(force(window),'b')
% plot(forcef(window),'r','linewidth',3)
% xlim([1 999])
% ylim([-.08 .08])
% xlabel('Time (ms)')
% ylabel('Measured Lift Force (N)')
% width = 5.8;height = 3;
% set(gca,'units','inches','position',[1 1 width height])
% set(gcf,'units','inches','position',[-14 1 width+2 height+1.5])
% set(gcf,'papersize',[width+1.5 height+1],'paperpositionmode','auto')
% saveFolder = '\\zeus\shared\Charlie Richter\Thesis\figures\';
% % saveas(gcf,[saveFolder,'force1'],'pdf')



%% Separate Data Ranges
regData = [];
for p = 5%1:size(fwin,1)
    domain = fwin(p,1):fwin(p,2);    % domain of interest
    [~,klocs1] = findpeaks(max(-kinematics(domain,1),0),...
        'minpeakheight',0,'minpeakdistance',100);%kpks = -kpks; % peaks in stroke angle
    klocs = klocs1+min(domain)-1;       % shift peak locations to true indices

    % Obtain flapping frequency in domain of interest
    Fs = 1000;L = length(kinematics(domain,1)); % sampling freq. and signal length
    NFFT = 2^nextpow2(L); % Next power of 2 from length of y
    Yfft = fft(kinematics(domain,1),NFFT)/L; % Fast Fourier Transform
    psd = 2*abs(Yfft(1:NFFT/2+1)); % Power Spectrum of Signal
    f = Fs/2*linspace(0,1,NFFT/2+1); % Associated frequency vector
    [~,psdLocs] = sort(findpeaks(psd),'descend');%,'sortstr','descend');
    flapfreq = f(psdLocs(1));
        
    amp = range(kinematics(domain,1));
    Re = 2*amp*flapfreq*spanTotal*chordTotal/15;
       
    % Lowpass filter force and kinematic data
##    cutoff = flapfreq*kBandwidth;width = 4;
##    h=fdesign.lowpass('Fp,Fst,Ap,Ast',cutoff,cutoff+width,.01,60,1000);
##    d=design(h,'equiripple'); %Lowpass FIR filter
##    forcef = filtfilt(d.Numerator,1,force); %zero-phase filtering
##    kf = filtfilt(d.Numerator,1,kinematics); %zero-phase filtering
    
    %% Octave does not have fdesign, so we need to use butter.
    [b, a]=butter(12, .048); % Best hand-tuned approximation of MATLAB design above.
    forcef = filtfilt(b,a,force);
    kf = filtfilt(b,a,kinematics);
    
    % Compute PSD of flapping in region
    Yforce = fft(force(domain),NFFT)/L;psdforce = 2*abs(Yforce(1:NFFT/2+1));
    Yforcef = fft(forcef(domain),NFFT)/L;psdforcef = 2*abs(Yforcef(1:NFFT/2+1));
    
    % Register data
    nCycles = length(klocs)-1;
    cycleLength1 = round(mean(diff(klocs)));
    tailLength = round(cycleLength1/2);
    cycleLength = cycleLength1+tailLength*2;
    meanDataHolder = nan(cycleLength,13,nCycles);

    whichvar = 11;
    for q = 1:nCycles;
        subwindow = (klocs(q)-tailLength:klocs(q)-tailLength+cycleLength-1);
        meanDataHolder(:,:,q) = [forcef(subwindow) kf(subwindow,:)]; % CHANGED 'KINEMATICS' TO 'KF' HERE
        h1=figure(1);hold on
        plot(rad2deg(meanDataHolder(:,whichvar,q)),'b--','linewidth',2)% 1 for force, 2 for stroke, 5 for def, 11 for dev
        drawnow
    end
    meanData = mean(meanDataHolder,3);
    dataLine = [meanData repmat([flapfreq Re kBandwidth],size(meanData,1),1)];
    regData = [regData;dataLine;nan(1,size(dataLine,2))]; %#ok<AGROW> 

    strokedot = meanData(:,3);
    shadecolor = .9*[1,1,1];
    
    axbuf = 1.1;
    ylimit = axbuf*[rad2deg(min(meanData(:,whichvar))) rad2deg(max(meanData(:,whichvar)))];
    cross = find([0;diff(sign(strokedot))]);
    fill([cross(1) cross(2)-1 cross(2)-1 cross(1)],100*[-1 -1 1 1],shadecolor,'edgecolor',shadecolor)
    fill([cross(3) cross(4)-1 cross(4)-1 cross(3)],100*[-1 -1 1 1],shadecolor,'edgecolor',shadecolor)
    
    for q = 1:nCycles;
        subwindow = (klocs(q)-tailLength:klocs(q)-tailLength+cycleLength-1);
        meanDataHolder(:,:,q) = [forcef(subwindow) kf(subwindow,:)]; % CHANGED 'KINEMATICS' TO 'KF' HERE
        h1=figure(1);hold on
        plot(rad2deg(meanDataHolder(:,whichvar,q)),'b--','linewidth',2)% 1 for force, 2 for stroke, 5 for def, 11 for dev
        drawnow
    end
    
    plot(rad2deg(meanData(:,whichvar)),'r','linewidth',1.5)

    xlabel('Time (ms)')
    ylim(ylimit)
    xlim([cross(1) 644])
    width = 5.8;height = 2;
    set(gcf,'units','inches','papersize',[width height],'position',[-8 6 width height],'paperpositionmode','auto')
    set(gca,'box','on','layer','top')
%     set(gca,'units','inches','position',[1 1 width height])
%     set(h1,'units','inches','position',[-14 1 width+2 height+2])
%     set(h1,'papersize',[width+1.5 height+1],'paperpositionmode','auto')
##    saveFolder = '\\zeus\shared\Charlie Richter\Thesis\figures\';
    saveFolder = '/home/crichter/ccsl/NSR20_Richter/paper/figures/';

%     ylabel('Measured Lift Force (N)')
%     ylabel('Stroke Angle \phi (\circ)')
%     ylabel('Deflection Angle \psi (\circ)')
    ylabel('Deviation Angle \theta (\circ)')
    
%     saveas(h1,[saveFolder,'regforce1'],'pdf')
%     saveas(h1,[saveFolder,'regstroke1'],'pdf')
%     saveas(h1,[saveFolder,'regdef1'],'pdf')
    saveas(h1,[saveFolder,'regdev1'],'pdf')

    %%
end

%% Export Data
% text = {'force','stroke','strokedot','strokeddot','def','defdot','defddot','cm','cmdot','cmddot','dev','devdot','devddot','freq','Re','kband','span','innerspan','totalspan','lowerchord','upperchord','totalchord','Xcm','Ycm','Zcm','Ixx','Ixy','Ixz','Iyx','Iyy','Iyz','Izx','Izy','Izz'};
% xlswrite([folderName 'eqdata.xls'],[text;text],1,'A1');
% geopropsarray = repmat([span innerSpan spanTotal...
%     chordLower chordUpper chordTotal...
%     X Y Z Ixx Ixy Ixz Iyx Iyy Iyz Izx Izy Izz],size(regData,1),1);
% geopropsarray(isnan(regData(:,1)),:) = nan;
% xlswrite([folderName 'eqdata.xls'],[regData geopropsarray],1,'A3');

% if toggleEndPlots
%     figure(),hold on
%     plot(100*regData(:,1),'b','linewidth',2)
%     plot(-1e-4*regData(:,10),'r','linewidth',2)
% end
% end