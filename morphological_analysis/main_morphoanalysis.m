%% Main script for FECG morphological analysis
% 
% This is the mains script for testing the morphological consistency of
% extracting the foetal signal using various methods. 
% Used extraction methods:
%  - ICA
%  - PCA
%  -piCA
% 
% Used morphological measures:
% - T/QRS ratio
% - ST segment
% - QT interval
% 
% 
% 
% NI-FECG simulator toolbox, version 1.0, February 2014
% Released under the GNU General Public License
%
% Copyright (C) 2014  Joachim Behar & Fernando Andreotti
% Oxford university, Intelligent Patient Monitoring Group - Oxford 2014
% joachim.behar@eng.ox.ac.uk, fernando.andreotti@mailbox.tu-dresden.de
%
% Last updated : 03-06-2014
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

%% Input parameters
% saving path
if isunix
%     path = ['/media/fernando/Data/foobar/fecgdata_test/' datestr(date,'yyyy.mm.dd') '/'];
 path = '/media/fernando/FetalEKG/2014.06_fecgsyn_simulations';
else
    path = ['C:\foobar\fecgdata\' datestr(date,'yyyy.mm.dd') '\'];
end

%% Set-up parameters
generate = 0;   % boolean, data should be generated? 
                % If not, path should direct to data location

% channels to be used in ICA                
% ch = 1:32;      
ch = [1:2:8 10:2:16 17:2:24 26:2:32];
debug = 0;

%% Data Generation
if generate
    mkdir(path)
    generate_data(path)  % generates set of unique simulated data for testing
else
    cd(path)
end

%% Extraction Methods
fls = dir('*.mat');     % looking for .mat (creating index)
fls =  arrayfun(@(x)x.name,fls,'UniformOutput',false);
stats_ica = zeros(length(fls),4);
stats_ts = zeros(length(fls),4);
tic
for i = 1:600   
    disp(['Extracting file ' fls{i} '..'])
    % = loading data
    load(fls{i})
    disp(num2str(i))
    noise = sum(cat(3,out.noise{:}),3);
    if isempty(noise)
        noise = zeros(size(out.mecg));
    end
    fs = out.param.fs;
    INTERV = round(0.05*fs);    % BxB acceptance interval
    TH = 0.3;                   % detector threshold
    REFRAC = round(.15*fs)/1000; % detector refractory period
    mixture = double(out.mecg) + sum(cat(3,out.fecg{:}),3) ...
        + noise;     % re-creating abdominal mixture
       
    
    %% Experiment 1
    % = preprocessing channels
    HF_CUT = 100; % high cut frequency
    LF_CUT = 0.7; % low cut frequency
    wo = 60/(fs/2); bw = wo/35;
    [b_lp,a_lp] = butter(5,HF_CUT/(fs/2),'low');
    [b_bas,a_bas] = butter(3,LF_CUT/(fs/2),'high');
    for j=ch
        lpmix = filtfilt(b_lp,a_lp,mixture(j,:));
        mixture(j,:) = filtfilt(b_bas,a_bas,lpmix);
    end
    % normalizing (significant for ICA)
    mixture = diag(1./max(mixture'))*mixture;
    
    % = using ICA
    disp('ICA extraction ..')
    loopsec = 60;   % in seconds
    icasig = ica_extraction(mixture,fs,ch,out.fqrs{1},loopsec);     % extract using IC
    
    % Calculate quality measures
    qrsica = qrs_detect(icasig,TH,REFRAC,fs);
    if isempty(qrsica)
        F1= 0;
        RMS = NaN;
        PPV = 0;
        SE = 0;
    else
        [F1,RMS,PPV,SE] = Bxb_compare(out.fqrs{1},qrsica,INTERV);
    end
    stats_ica(i,:) = [F1,RMS,PPV,SE];
    
    % = using TSc
    disp('TS extraction ..')
    % look for channel with largest SNRfm
    amps = sum(double(out.fecg{1}).^2,2)./sum(double(out.mecg).^2,2);
    [~,chts]=max(amps);       % chosing the channel with highest fetal signal ratio
    residual = mecg_cancellation(out.mqrs,mixture(chts,:),'TS-CERUTTI',0);
    qrsts = qrs_detect(residual,TH,REFRAC,fs);
    [F1,RMS,PPV,SE] = Bxb_compare(out.fqrs{1},qrsts,INTERV);
    stats_ts(i,:) = [F1,RMS,PPV,SE];

    % Debug plots
    if debug && ~isempty(qrsica)
        close all
        FONT_SIZE = 15;
        LINE = 2;
        MSIZE = 7;
        tm = 1/fs:1/fs:length(residual)/fs;
        figure('name','MECG cancellation');
        subplot(1,2,1)
        plot(tm,mixture(chts,:),'k','LineWidth',LINE-1);
        hold on, plot(tm,icasig-2,'-b','LineWidth',LINE);
        plot(out.fqrs{1}/fs,-1,'xr','MarkerSize',MSIZE)
        plot(tm(qrsica),-1,'or','LineWidth',LINE,'MarkerSize',MSIZE);
        legend('mixture','residual','reference FQRS','detected FQRS');
        title('ICA for extracting the FECG');
        xlabel('Time [sec]'); ylabel('Amplitude [NU]')
        
        subplot(1,2,2)
        
        plot(tm,mixture(chts,:),'k','LineWidth',LINE-1);
        hold on, plot(tm,residual-2,'b','LineWidth',LINE);
        plot(out.fqrs{1}/fs,-1,'xr','MarkerSize',MSIZE)
        plot(tm(qrsts),-1,'or','LineWidth',LINE,'MarkerSize',MSIZE);
        legend('mixture','residual','reference FQRS','detected FQRS');
        title('Template subtraction for extracting the FECG');
        xlabel('Time [sec]'); ylabel('Amplitude [NU]')
        set(findall(gcf,'type','text'),'fontSize',FONT_SIZE);
        linkaxes


    end
    clearvars -except stats_ica stats_ts fls ch debug

end
toc
%% Morphological Analysis


