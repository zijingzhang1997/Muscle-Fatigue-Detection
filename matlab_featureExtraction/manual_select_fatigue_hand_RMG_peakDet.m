%% calculate max correlation time lag, 
% DTW another option, not always as good as time lag (s), no output of time
% lag

close all
clear all
Casenum='2';    

RoutineNList={'grasp','p23','gripper','fingerInst','GripRing'...
    };
RoutineN=RoutineNList{1};
% RoutineN='grip';
RNum=2;
saveOPt=1;
ReloadOpt=0;

name=[Casenum,RoutineN];
Tname=strrep(name,'_',' ');
ExpDate='fatigue_ZZ_hand_new_v1';
%ExpDate='memory_ges_v1_ZZ';
ExpDate='fatigue_TC_hand_v1';
dataPath=['D:\legRMG\data\',ExpDate,'\'];
addpath('D:\eye RMG\matlab code_all');
addpath('D:\COVID\COVID_HF_spectrum\dyspnea_study_new_code')

fs=5e3;
fsDS=500;
CaseName=['Case',name];
fileName=[CaseName,'Routine',num2str(RNum)];
filt=[0.05,8];
% make sure filter LF is lower so higher frequency is filtered out

% Less noisy for peak detection

%toff=[10:length(Ch_data_raw)/fs-2]';

filePathName = [dataPath,fileName,'.tdms'];
filePathName_m = [dataPath,fileName,'.mat'];
if ~exist(filePathName_m,'file')
   convertTDMS(true,filePathName);
end
load(filePathName_m);


load([dataPath,fileName,'.mat']);

opt.filtType = 'LpHp'; opt.orderHP = 5;
opt.f3db = filt(1); opt.fpLP = filt(2); opt.fstLP = opt.fpLP+1;
Ch_num=cat(1,[3:18]',[23:38]');




for i = 1:32
  
  Ch_data(:,i)=ConvertedData.Data.MeasuredData(Ch_num(i)).Data; % ch1 amp
  
  Ch_data_raw=Ch_data;
    if mod(i,2) == 0  % if phase channel 
      Ch_data(:,i)=unwrap(deg2rad(Ch_data(:,i)));
    end
toff=[10:length(Ch_data_raw)/fs-15]';  
% toff=[30:length(Ch_data_raw)/fs-20]'; 
% % %toff=[5:55]';  
%toff=[10+180:170+180]';  
  [Ch_data_filt(:,i),Ch_data_filt_raw(:,i)]= filtSignal(Ch_data(:,i),toff,opt,fs,fsDS);
 % Ch_data_filt=Ch_data_filt_raw;
  Ch_data_filt(:,i)=detrend(Ch_data_filt(:,i));
end 

%%  EOG 
% EOG is less noisy, can directly use the same filter of low pass 
% so noisy peaks can be smoothed , better for peak detection 
filt=[0.1,8];
opt_emg.filtType = 'LpHp'; opt_emg.orderHP = 5;
opt_emg.f3db = filt(1); opt_emg.fpLP = filt(2); opt_emg.fstLP =opt_emg.fpLP+1;
opt_emg.SMorder=5;opt_emg.SMlen=9;  %% smoothing filter parameter , order less, length longer, more smooth 
opt_emg.movmean=50;
opt_emg.fsDS=fsDS;
% opt_emg.SMlen=199; opt_emg.movmean=100;
Ch_num_emg=[19,20];
for i = 1:2  %% use first two channels of EMG 
  
  Ch_data_emg(:,i)=ConvertedData.Data.MeasuredData(Ch_num_emg(i)).Data; % ch1 amp
  
  [Ch_data_emg_filt(:,i),Ch_data_emg_filt_raw(:,i),Ch_data_emg_en(:,i),Ch_data_emg_raw(:,i)]= filtSignal_emg(Ch_data_emg(:,i),toff,opt_emg);
  % Ch_data_emg_filt_raw : only filter, no smooth 
  % Ch_data_emg_filt  : filter + smooth + normalization
%  Ch_data_emg_filt=Ch_data_emg_filt_raw;
  Ch_data_emg_filt(:,i)=detrend(Ch_data_emg_filt(:,i));
  
end 


Ch_data_filt=[Ch_data_filt,Ch_data_emg_filt];
Ch_data_filt = sgolayfilt(Ch_data_filt,3,71);
Ch_data_emg=cat(3,Ch_data_emg_raw,Ch_data_emg_en,Ch_data_emg_filt_raw,Ch_data_emg_filt);

Ch_data_filt=detrend(Ch_data_filt);
Chan_Name={'Tx1Rx1 amp','TxRx1 ph','Tx2Rx1 amp','Tx2Rx1 ph','Tx3Rx1 amp','Tx3Rx1 ph','Tx4Rx1 amp','Tx4Rx1 ph',...
    'Tx1Rx2 amp','Tx1Rx2 ph','Tx2Rx2 amp','Tx2Rx2 ph','Tx3Rx2 amp','Tx3Rx2 ph','Tx4Rx2 amp','Tx4Rx2 ph',...
    'Tx1Rx3 amp','TxRx3 ph','Tx2Rx3 amp','Tx2Rx3 ph','Tx3Rx3 amp','Tx3Rx3 ph','Tx4Rx3 amp','Tx4Rx3 ph',...
    'Tx1Rx4 amp','TxRx4 ph','Tx2Rx4 amp','Tx2Rx4 ph','Tx3Rx4 amp','Tx3Rx4 ph','Tx4Rx4 amp','Tx4Rx4 ph','EMG1','EMG2'}';
ind=[1,2,11,12,21,22,31,32,33,34];
ind2=[3,4,9,10,23,24,29,30,33,34];
% 8 slef channels , change direction - + 
% because BR is calculated based on  Max peak

%% frequency calculation
t=1:round(max(toff)-min(toff));
MotionTrue=zeros(length(t),1);
MotionTrue(:)=60;
Motion=[t'  MotionTrue];
% 8 slef channels , change direction - + 
% because BR is calculated based on  Max peak
opts3.tWinBR = 4; % Window on which br is estimated

opts3.minInterceptDist = 0.2; % minimum time (s) between two intercepts default 0.05
opts3.calibPk = 1; % Calibrate out peaks of small height (max-min)
opts3.calibMinPkRatio = 0.3;
opts3.tWinVar=10;opts3.tWinpp=10;
opts3.tWinPk=5; % Window for peak detection moving average
% select the best channel by min cov
ind_rmg=[ind(1:8) ind2(1:8)];
opts3.fig=0;


h_w1=plotWaveform_s(Ch_data_filt(:,ind),Chan_Name(ind),ind,fsDS,[Tname,'self CH']);

h_w2=plotWaveform_s(Ch_data_filt(:,ind2),Chan_Name(ind2),ind2,fsDS,[Tname,'cross CH']);



%%feat
% change channel manually find the best channel 


 sign_opt=-1;
RMG_ch_opt=9;
EMG_ch_opt=33;
RMG_ch_featDet1=24;
RMG_ch_featDet2=21;
% if reload channel number 
if ReloadOpt==1
load([dataPath,'mat','\',fileName,'\VarFeat','.mat'],'RMG_ch_opt','EMG_ch_opt',...
  'toff','sign_opt','RMG_ch_featDet1','RMG_ch_featDet1');
fprintf('reload opt ch');
end

% sign_opt is the RMG sign + - for the optimal channel selected, will be
% reloaded and remebered 

% if override reloadded channel 
% sign_opt=1;
% RMG_ch_opt=21;
% EMG_ch_opt=34;
Ch_data_filt(:,RMG_ch_opt)=Ch_data_filt(:,RMG_ch_opt)* sign_opt;
%% for feature detection assign two RMG channels (not neccessarily the same as Opt for corelation)

opts3.fig=1;
% feat ; BR PP IN EX (window mean) + BR PP IN EX (window var)  size same as input 
% featMean BR PP IN EX mean of (window mean) + BR PP IN EX mean of (window var) 
[feat_rmg1,pk_rmg1,featMean_rmg1,h_rmg(1)] = brEstAvg(Ch_data_filt(:,RMG_ch_featDet1),fsDS,opts3,Chan_Name(RMG_ch_featDet1),Motion,100);
[feat_rmg2,pk_rmg2,featMean_rmg2,h_rmg(2)] = brEstAvg(Ch_data_filt(:,RMG_ch_featDet2),fsDS,opts3,Chan_Name(RMG_ch_featDet2),Motion,100);
[feat_emgOpt,pk_emgOpt,featMean_emgOpt,h_emgOpt] = brEstAvg(Ch_data_filt(:,EMG_ch_opt),fsDS,opts3,Chan_Name(EMG_ch_opt),Motion,100);



%% EMG RMG correlation
% for all ime 

tOff=((0:(length(Ch_data_filt)-1))/fsDS)';
% OptFig=1;
% [max_coef,max_r,h_cor,h_overlap,dist]=max_cor(Ch_data_filt,[1/fsDS,max(tOff)],fsDS,RMG_ch_opt,EMG_ch_opt,OptFig,Chan_Name);
%% time window correlation 
tWin=5; tSlide=2;
% tWin=15; tSlide=5;
% tWin=18; tSlide=3;
tStart=1/fsDS;

% start / end time of each epoch  2* epochNum
tEnd=toff(end)-toff(1)-tStart-tWin; 
StartEndT=cat(1,tStart:tSlide:tEnd, tStart+tWin:tSlide:tEnd+tWin); 
max_coef_seg=[];
max_r_seg=[];
dist=[];
for i=1:size(StartEndT,2)
    [StartEndT(1,i),StartEndT(2,i)]
[max_coef_seg(i),max_r_seg(i),~,~]=max_cor(Ch_data_filt,[StartEndT(1,i),StartEndT(2,i)],fsDS,RMG_ch_opt,EMG_ch_opt,0,Chan_Name);



end


% remove points with tLag< 0  
% linear fitting 

Correct_T=StartEndT(1,(max_r_seg>0));
Correct_max_r_seg=max_r_seg(max_r_seg>0)*1e3;
[p,S] = polyfit(Correct_T,Correct_max_r_seg,1); 
[y_fit,delta] = polyval(p,Correct_T,S);




sz=13;
h_segCor=figure()
scatter(StartEndT(1,:),max_r_seg,'o')
xlabel('Epoch Start Time(s)','FontSize',sz)
ylabel('EMG&RMG Time Lag(s)','FontSize',sz)
% ylim([0 max(max_r_seg)+0.02])
title([Tname,' ',Chan_Name{RMG_ch_opt},' ',Chan_Name{EMG_ch_opt},...
    '  std=', num2str(std(Correct_max_r_seg),2),'  mean=', num2str(mean(Correct_max_r_seg),2)])

h_segCor_fit=figure()
scatter(Correct_T,Correct_max_r_seg,'o')
xlabel('Epoch Start Time(s)','FontSize',sz)
ylabel('EMG&RMG Time Lag(ms)','FontSize',sz)
%ylim([0 max(Correct_max_r_seg)+20])
hold on 
plot(Correct_T,y_fit,'r--')
plot(Correct_T,y_fit+2*delta,'m--',Correct_T,y_fit-2*delta,'m--')
title('Linear Fit of Data with 95% Prediction Interval')
legend('Data','Linear Fit','95% Prediction Interval','Location','southeast')
txtData=[' \mu T_{lag} =', num2str(mean(Correct_max_r_seg),3),' ms;',' \sigma T_{lag} =', num2str(std(Correct_max_r_seg),3),' ms;',...
    '\Delta mean=',num2str(mean(delta),3),' ms;',' p_{1}=',num2str(p(1),2)];
title({[Tname,' ',Chan_Name{RMG_ch_opt},' ',Chan_Name{EMG_ch_opt}],txtData},'FontSize',sz);

%% save figures 

status = mkdir([dataPath,'mat\',fileName,'\']);
status = mkdir([dataPath,'fig\',fileName,'\']);
saveFigFolder=[dataPath,'fig\',fileName,'\'];
saveMatFolder=[dataPath,'mat\',fileName,'\'];

if saveOPt==1
  
figName = [saveFigFolder,'waveform'];
print(h_w1,[figName,'.tiff'],'-dtiff','-r300');
savefig(h_w1,[figName,'.fig']);

figName = [saveFigFolder,'waveform_crossCH'];
print(h_w2,[figName,'.tiff'],'-dtiff','-r300');
savefig(h_w2,[figName,'.fig']);
    
    
figName = [saveFigFolder,'rmg','_ch_',num2str(RMG_ch_featDet1)];
print(h_rmg(1),[figName,'.tiff'],'-dtiff','-r300');
savefig(h_rmg(1),[figName,'.fig']);
figName = [saveFigFolder,'rmg','_ch_',num2str(RMG_ch_featDet2)];
print(h_rmg(2),[figName,'.tiff'],'-dtiff','-r300');
savefig(h_rmg(2),[figName,'.fig']);
figName = [saveFigFolder,'emg','_ch_',num2str(EMG_ch_opt)];
print(h_emgOpt,[figName,'.tiff'],'-dtiff','-r300');
savefig(h_emgOpt,[figName,'.fig']);   
        
    
% figName = [saveFigFolder,'_Cor_Fig','_ch',num2str(RMG_ch_opt),'_',num2str(EMG_ch_opt)];
% print(h_cor,[figName,'.tiff'],'-dtiff','-r300');
% savefig(h_cor,[figName,'.fig']);
% figName = [saveFigFolder,'overlap_Fig','_ch',num2str(RMG_ch_opt),'_',num2str(EMG_ch_opt)];
% print(h_overlap,[figName,'.tiff'],'-dtiff','-r300');
% savefig(h_overlap,[figName,'.fig']);

figName = [saveFigFolder,'segCor_Fig','_ch',num2str(RMG_ch_opt),'_',num2str(EMG_ch_opt)];
print(h_segCor,[figName,'.tiff'],'-dtiff','-r300');
savefig(h_segCor,[figName,'.fig']);
figName = [saveFigFolder,'segCor_fit_Fig','_ch',num2str(RMG_ch_opt),'_',num2str(EMG_ch_opt)];
print(h_segCor_fit,[figName,'.tiff'],'-dtiff','-r300');
savefig(h_segCor_fit,[figName,'.fig']);


save([saveMatFolder,'VarFeat','.mat'],...
    'feat_rmg1','feat_rmg2','feat_emgOpt','RMG_ch_featDet1','RMG_ch_featDet2','RMG_ch_opt','EMG_ch_opt',...
    'max_coef_seg','max_r_seg','StartEndT','toff','sign_opt','Correct_T','p','delta',...
    'Correct_max_r_seg','y_fit');
end
%%

function [ampCh_filt_norm,ampCh_filt]=filtSignal(ampCh,toff,opt,fs,fsDS)
    
   
  
    ampCh=resample(ampCh,fsDS,fs);
    ampCh=ampCh((toff(1)*fsDS):toff(size(toff))*fsDS);
    ampCh_filt = filterLpHp(ampCh,fsDS,opt); % th amp
    
    ampCh_filt_norm = normalize(ampCh_filt);

end 

function [ampCh_filt_norm,ampCh_filt,ampCh_en,ampCh_raw]=filtSignal_emg(ampCh,toff,opt)
    
    fs=1e3;  % biopac sampling rate =1e3 
    fsDS=opt.fsDS; % down sampling rate  same as NCS  
    
  
    ampCh=resample(ampCh,fsDS,fs);
    ampCh_raw=ampCh((toff(1)*fsDS):toff(size(toff))*fsDS);
    
    [ampCh_en,~] = envelope(ampCh_raw,5,'peak');
    
    ampCh_en=detrend(ampCh_en);
    ampCh_filt = filterLpHp(ampCh_en,fsDS,opt); % th amp
    
    
    ampCh_filt_sm = sgolayfilt(ampCh_filt,opt.SMorder,opt.SMlen);
    %ampCh_filt_sm = smoothdata(ampCh_filt,'movmean',200);
    
    ampCh_filt_norm = normalize(ampCh_filt_sm);
    
end 

