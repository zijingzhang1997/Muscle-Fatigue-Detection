%revise by zijing from 'brEst'
% BR is average of all BRs over the time window, not pick only one BR
%calculate variation of BR, PP, inhale /exhale time
function [feat,pk,feature_all,h] = brEstAvg(Data,fs,opts,Chan_Name,Motion,fsDS2)

if ~isfield(opts,'tWinBR')
    opts.tWinBR = 8;
 %   fprintf('Default BR estimation window: %3.2f\n',opts.tWinBR);
end
if ~isfield(opts,'calibPk')
    opts.calibPk = 0; % Do not calibrate by default
else
    if ~isfield(opts,'calibT')
        opts.calibT = [20,40]; 
      %  fprintf('Default BR estimation peak calibration window [%d, %d].\n',opts.calibT(1),opts.calibT(2));
    end
    if ~isfield(opts,'calibMinPkRatio')
        opts.calibMinPkRatio = 0.4;
       % fprintf('Default min peak height can be %d%% of avg peak height in calibration window.\n',100*opts.calibMinPkRatio);
    end
end
Data=resample(Data,fsDS2,fs);
fs=fsDS2;
t = (0:(length(Data)-1))/fs;

% -------------------------------------------------------------------------
% Minima and maxima detection on NCS thorax and abdomen data
% -------------------------------------------------------------------------
pk = findMaxMin(Data,fs,opts);

if pk(1).ind(1) == 1
    % If first peak is inhalation peak, skip it
    pk(1).ind = pk(1).ind(2:end); 
    pk(1).idx = pk(1).idx(2:end);
end

if pk(1).ind(end) == 1
    pk(1).ind = pk(1).ind(1:end-1);
    pk(1).idx = pk(1).idx(1:end-1);
end


pkMax = pk(1).idx(pk(1).ind == 1);
pkMin = pk(1).idx(pk(1).ind == 0);

if opts.calibPk == 1
  %  fprintf('\nBR: Performing calibration of pk-pk height.\n')
    % This is the change in Pk-Pk thorax and abdomen signal
    del2 = zeros(length(pkMax),1);

    % So for a cycle, considering there exists 2 minima and 2 maxima point:
    % Calculation is peformed using the difference between maxima and first
    % minima. The update is performed at the end of cycle to be consistent with
    % TV from airflow calculation. And this value is held until next update or
    % the end of the waveform.
    for i = 1:length(pkMax)
            del2(i) = abs(Data(pkMax(i),1)-Data(pkMin(i),1)); % Making it positive always
    end

    tBRmax1 = t(pkMax);
    idx2Calib = ((tBRmax1 >= opts.calibT(1))&(tBRmax1 <= opts.calibT(2)));
    delNcs2Calib = mean(del2(idx2Calib)); 
%     fprintf('delNCS2Calib = %f\n',delNcs2Calib);

    idxPk = (del2 >= opts.calibMinPkRatio*delNcs2Calib);

    pkMax = pkMax(idxPk); % Only keep indices that satisfy, otherwise ignore
    pkMin = pkMin(idxPk);
    
    pk(1).idxValidPk = idxPk;
end

% fig1 = figure;  by zijing
% nFig = size(ncsData,2)+1;

%%

%% Only one NCS column    
    % Find breath rate in last tWinBR sec (or slightly less)
    mean1=zeros(length(t),1);
    mean2=zeros(length(t),1);
    mean3=zeros(length(t),1);
    mean4=zeros(length(t),1);
    br = zeros(length(t),1);
    var1=zeros(length(t),1);  % BR var
    var2=zeros(length(t),1);  %pp var
    var3=zeros(length(t),2);  %inhale exhale time var
  
    t = t(:);
    tBRmax1 = t(pkMax);

    for iter = 1:length(t)

        idxBR = find((tBRmax1 > (t(iter)-opts.tWinBR))&(tBRmax1 <= t(iter)));

        if (size(idxBR,1) <= 1)
            br(iter,1) = 0;
            
        else  
            %average calculatio of BR
            for i=1:length(idxBR)-1
               br(iter,1) = br(iter,1)+(idxBR(i+1)- idxBR(i))/(tBRmax1(idxBR(i+1))-tBRmax1(idxBR(i)))/(length(idxBR)-1);   
            end       
            % br(iter,1) = (idxBR(2)-idxBR(1))/(tBRmax1(idxBR(2))-tBRmax1(idxBR(1)));
            %only use one cycle    
        end
    end
    %BR variation
    for iter = 1:length(t)

        idxBR = find((tBRmax1 > (t(iter)-opts.tWinVar))&(tBRmax1 <= t(iter)));

        if (size(idxBR,1) <= 1)
            var1(iter,1) = 0;
        else  
            %average calculatio of BR
            temp1=zeros(length(idxBR)-1,1);
            
            for i=1:length(idxBR)-1
               temp1(i) = (idxBR(i+1)- idxBR(i))/(tBRmax1(idxBR(i+1))-tBRmax1(idxBR(i))); 
               temp1(i)=temp1(i)*60;
            end
            var1(iter,1)=var(temp1)/((mean(temp1)).^2);
            mean1(iter,1)=mean(temp1);
            % br(iter,1) = (idxBR(2)-idxBR(1))/(tBRmax1(idxBR(2))-tBRmax1(idxBR(1)));
            %only use one cycle    
        end
    end 
    %peak-peak variation
       for iter = 1:length(t)

        idxBR = find((tBRmax1 > (t(iter)-opts.tWinpp))&(tBRmax1 <= t(iter)));

        if (size(idxBR,1) <= 1)
            var2(iter,1) = 0;
            
            
        else  
            %average calculatio of BR
            temp1=zeros(2*(length(idxBR)-1),1);
            temp_in=zeros(length(idxBR)-1,1);
            temp_ex=zeros(length(idxBR)-1,1);
            for i=1:length(idxBR)-1
                a1=pkMax(idxBR(i));
                a2=pkMax(idxBR(i+1));
                b=pkMin(idxBR(i+1));
                temp1(2*i-1) = Data(a2)-Data(b);
                temp_in(i)=t(a2)-t(b);
                temp_ex(i)=t(b)-t(a1);
                temp1(2*i) = Data(a1)-Data(b); 
               
            end
            mean2(iter,1)=mean(temp1);mean3(iter,1)=mean(temp_in);mean4(iter,1)=mean(temp_ex);
            var2(iter,1)=var(temp1)/((mean(temp1))^2);%CoV = std/mean
            var3(iter,1)=var(temp_in)/((mean3(iter,1))^2);  %normalize??
            var3(iter,2)=var(temp_ex)/((mean4(iter,1))^2);
        end
        end
      
 var2=var2.*100; 
 var1=var1.*100; 
 var3=var3.*100; 
 Var=[var1,var2,var3];
    br = 60.*br;
    
meanBR=mean(mean1);meanPP=mean(mean2);meanIn=mean(mean3);meanEx=mean(mean4);
VarBR_mean=mean(var1(var1~=0));
VarBR_std=std(var1(var1~=0));
VarPP_mean=mean(var2(var2~=0));
VarPP_std=std(var2(var2~=0));
VarIn=var3(:,1);
VarIn_mean=mean(VarIn(VarIn~=0));
VarIn_std=std(VarIn(VarIn~=0));
VarEx=var3(:,2);
VarEx_mean=mean(VarEx(VarEx~=0));
VarEx_std=std(VarEx(VarEx~=0));

VarFeature=[ meanBR; meanPP; meanIn; meanEx;VarBR_mean; VarPP_mean; VarIn_mean; VarEx_mean];
feat=[mean1 mean2 mean3 mean4 var1 var2 var3(:,1) var3(:,2)];
%%
    %% Find features in each cycle (min-max-min)
cycle=length(pkMax)-1;
br = zeros(cycle,1);ibi = zeros(cycle,1);pp = zeros(cycle,1);in = zeros(cycle,1);ex = zeros(cycle,1);
ibi = zeros(cycle,1);  % inter breath time interval 
IEpp= zeros(cycle,1);  % inhale volume / exhale volume 
IER= zeros(cycle,1);  % inhale time / exhale time ratioe 
skew=zeros(cycle,1);kurt=zeros(cycle,1);
for i=1:cycle
        br(i)=60/(t(pkMin(i+1))-t(pkMin(i)));  %BR unit BPM 
        ibi(i)=t(pkMin(i+1))-t(pkMin(i));  %Inter-breath interval
        pp(i)=(Data(pkMax(i))-Data(pkMin(i))-Data(pkMin(i+1)))/2;  %% peak to peak averaged on min-max & max-min 
        IEpp(i) =(Data(pkMax(i))-Data(pkMin(i)))/(Data(pkMax(i))-Data(pkMin(i+1)));
    
        in(i)=t(pkMax(i))-t(pkMin(i));
        ex(i)=t(pkMin(i+1))-t(pkMax(i));
        IER(i)=in(i)/ex(i);
        
        
        DataSeg=Data(pkMin(i):pkMin(i+1));
        kurt(i)=kurtosis(DataSeg);  %measure of the "tailedness"
        skew(i)=skewness(DataSeg);  %measure of the asymmetry
        en(i)=entropy(DataSeg);
    end
 
entro=entropy(Data); 
skew_mean=mean(skew);kurt_mean=mean(kurt);

covBR=std(br)/mean(br);
covIBI=std(ibi)/mean(ibi);
covPP=std(pp)/mean(pp);
covIN=std(in)/mean(in);
covEX=std(ex)/mean(ex);

% auto correlation and successive difference in neighbor cycles 
Cor_br = xcorr(br,1,'coeff');Cor_br =Cor_br (1);
SD_br= mean(abs(diff(br))./br(1:end-1));
Cor_ibi = xcorr(ibi,1,'coeff');Cor_ibi =Cor_ibi (1);
SD_ibi= mean(abs(diff(ibi))./ibi(1:end-1));
Cor_pp = xcorr(pp,1,'coeff');Cor_pp =Cor_pp (1);
SD_pp= mean(abs(diff(pp))./pp(1:end-1));
Cor_in = xcorr(in,1,'coeff');Cor_in =Cor_in (1);
SD_in= mean(abs(diff(in))./in(1:end-1));
Cor_ex = xcorr(ex,1,'coeff');Cor_ex =Cor_ex (1);
SD_ex= mean(abs(diff(ex))./ex(1:end-1));
Cor_ex = xcorr(ex,1,'coeff');Cor_ex =Cor_ex (1);
SD_ex= mean(abs(diff(ex))./ex(1:end-1));
Cor_IEpp = xcorr(IEpp,1,'coeff');Cor_IEpp =Cor_IEpp (1);
SD_IEpp= mean(abs(diff(IEpp))./IEpp(1:end-1));
Cor_IER = xcorr(IER,1,'coeff');Cor_IER =Cor_IER (1);
SD_IER= mean(abs(diff(IER))./IER(1:end-1));

feature_cycle=[mean(br),std(br),mean(ibi),std(ibi),mean(pp),std(pp),mean(in),std(in),mean(ex),std(ex),...
    mean(IEpp),std(IEpp),mean(IER),std(IER),...
        covBR,covPP,covIN,covEX,covIBI,... 
        Cor_br, SD_br,Cor_ibi,SD_ibi,Cor_pp,SD_pp,Cor_in,SD_in,Cor_ex,SD_ex,Cor_IEpp,SD_IEpp,Cor_IER,SD_IER,...
        skew_mean, kurt_mean,entro,cycle];


%%
[FeatAll]=Cal_feat(feat,fs,opts,Motion);
VarFeature=[VarFeature; FeatAll'];

feature_all=[VarFeature;feature_cycle'];
feature_all(isnan(feature_all))=100;
%%
h=[];
if opts.fig==1
nFig=4;
    h=figure();
    ax(1) = subplot(nFig,1,1);
    plot(t,Data); hold on;
    plot(t(pkMax),Data(pkMax,1),'^',...
        t(pkMin),Data(pkMin,1),'v');
    leg = {'waveform','Max','Min'};
    
    
  plotCute1('Time (s)','amp(a.u.)',ax(1),[],leg,1,'Horizontal');
    
    title(Chan_Name);
    ax(2) = subplot(nFig,1,2);
    plot(t,mean1(:,1),'LineWidth',2);
    hold on 
    plot(Motion(:,1),Motion(:,2),'-.','LineWidth',2);
    leg = {'Freq','True'};
    plotCute1('Time (s)','Motion Freq (BPM)',ax(2),['mean(BPM):',num2str(meanBR,2),' Var(%):',num2str(VarBR_mean,2)],leg,1,'Horizontal');
    
     ax(3) = subplot(nFig,1,3);
    plot(t,mean2,'LineWidth',2);
    hold on 
    plotCute1('Time (s)','Mean PP',ax(3),[],[],[],'Horizontal');
         ax(4) = subplot(nFig,1,4);
    plot(t,var2,'LineWidth',2);
    hold on 
    plotCute1('Time (s)','Variation PP',ax(4),['mean(PP):',num2str(meanPP,2),' Var(%):',num2str(VarPP_mean,2)],[],[],'Horizontal');
    linkaxes(ax,'x');
    set(gcf,'Position',[200,100,1200,800]);
   
    
    
     end
end


function [FeatAll]=Cal_feat(Winfeat,fs,opt,Motion)

%feat=[mean1 mean2 mean3 mean4 var1 var2 var3(:,1) var3(:,2)]; 
% variation is in (%)
% featMean BR PP IN EX mean of (window mean) + BR PP IN EX mean of (window var) 
% feat ; BR PP IN EX (window mean) + BR PP IN EX (window var)  size same as input 

t = (0:(length(Winfeat)-1))/fs;




% truncate the beginning artificate use time tWinPk +0.5
tOff=opt.tWinPk+0.5;
t=t(tOff*fs:end);
Winfeat=Winfeat(tOff*fs:end,:);

featMean=mean(Winfeat,1);

% calculate linear fit of BR and pp slope *100
[p,S] = polyfit(t,Winfeat(:,2),1); 
[y_fit,delta] = polyval(p,t,S);
fit_parameter_pp=[p(1)*100,p(2),mean(delta)];
[p,S] = polyfit(t,Winfeat(:,1),1); 
[y_fit,delta] = polyval(p,t,S);
fit_parameter_br=[p(1)*100,p(2),mean(delta)];
% error of frequency  abs *100 
err_br_mean=mean(abs((Winfeat(:,1)-Motion(1,2))))/Motion(1,2)*100;


FeatAll=[fit_parameter_pp,fit_parameter_br,err_br_mean];



end
