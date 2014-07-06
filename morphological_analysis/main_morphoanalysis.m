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
 path = '/media/fernando/FetalEKG/2014.07_fecgsyn_simulations(2.0)/';
else
    path = ['C:\foobar\fecgdata\' datestr(date,'yyyy.mm.dd') '\'];
end

%% Set-up parameters
generate = 1;   % boolean, data should be generated? 
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
for i = 1:100   
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
    mixture = mixture(ch,:);    % reducing number of channels
    
    %% Experiment 1
    % = preprocessing channels
    HF_CUT = 100; % high cut frequency
    LF_CUT = 0.7; % low cut frequency
    wo = 60/(fs/2); bw = wo/35;
    [b_lp,a_lp] = butter(5,HF_CUT/(fs/2),'low');
    [b_bas,a_bas] = butter(3,LF_CUT/(fs/2),'high');
    for j=1:length(ch)
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
    for j = 1:length(ch)
        residual = mecg_cancellation(out.mqrs,mixture(j,:),'TS-CERUTTI',0);
        qrsts = qrs_detect(residual,TH,REFRAC,fs);
        [F1(j),RMS(j),PPV(j),SE(j)] = Bxb_compare(out.fqrs{1},qrsts,INTERV);
    end
    [~,maxch] = max(F1);
    stats_ts(i,:) = [F1(maxch),RMS(maxch),PPV(maxch),SE(maxch)];
    
    
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
if debug
    
    % Plot about cases
    
    % This script plots boxplots with 2 groups  
    % Boxplot multicolor
    stats_ica(:,[1 3 4]) = 100.*stats_ica(:,[1 3 4]);
    stats_ts(:,[1 3 4]) = 100.*stats_ts(:,[1 3 4]);
    LINE_WIDTH = 1.2;
    FONT_SIZE = 12;
    MARKER_SIZE = 7;
    fig1=figure(1)
    base = cellfun(@(x) ~isempty(regexp(x, '_l\d.mat$', 'match')), fls);   
    c1 = cellfun(@(x) ~isempty(regexp(x, '_c1.mat$', 'match')), fls);
    c2 = cellfun(@(x) ~isempty(regexp(x, '_c2.mat$', 'match')), fls);
    c3 = cellfun(@(x) ~isempty(regexp(x, '_c3.mat$', 'match')), fls);
    c4 = cellfun(@(x) ~isempty(regexp(x, '_c4.mat$', 'match')), fls);
    c5 = cellfun(@(x) ~isempty(regexp(x, '_c5.mat$', 'match')), fls);

    stats_ica = stats_ica(:,1);
    stats_ts = stats_ts(:,1);
    N = sum(base)+sum(c1)+sum(c2)+sum(c3)+sum(c4)+sum(c5);
    bh=boxplot([stats_ica(base); stats_ica(c1); stats_ica(c2); stats_ica(c3); stats_ica(c4); stats_ica(c5);...
        stats_ts(base); stats_ts(c1); stats_ts(c2); stats_ts(c3); stats_ts(c4); stats_ts(c5)], ...
        {[repmat({'Base'},1,sum(base)) repmat({'Case1'},1,sum(c1)) repmat({'Case2'},1,sum(c2)) repmat({'Case3'},1,sum(c3)) repmat({'Case4'},1,sum(c4)) repmat({'Case5'},1,sum(c5)) ...
        repmat({'Base'},1,sum(base)) repmat({'Case1'},1,sum(c1)) repmat({'Case2'},1,sum(c2)) repmat({'Case3'},1,sum(c3)) repmat({'Case4'},1,sum(c4)) repmat({'Case5'},1,sum(c5)) ] ...
        [repmat({'ICA'},N,1);repmat({'TS'},N,1)]},'colors',repmat('rk',1,6),'factorgap',10,'labelverbosity','minor','labelorientation','inline')
    set(bh,'linewidth',LINE_WIDTH);
    ylabel('F1 (in %)','FontSize',FONT_SIZE)
    xlabel('Recording','FontSize',FONT_SIZE)
    set(gca,'FontSize',FONT_SIZE)
    set(findall(gcf,'type','text'),'fontSize',FONT_SIZE);
%     ylim([0 105])
    xl=xlabel('Recording Number','FontSize',FONT_SIZE);
    pos=get(xl,'Pos');
    set(xl,'Pos',[pos(1) pos(2)-30 pos(3)])
    save2pdf('boxplot',fig1,600)
 
    % Plot about SNRmn
    a = cellfun(@(x) strsplit(x,'_snr'), fls,'UniformOutput',0);
    
    b = cellfun(@(x) length(x)>1,a);
end



%% Morphological Analysis


