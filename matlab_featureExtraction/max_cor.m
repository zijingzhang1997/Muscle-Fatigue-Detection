function [max_coef,max_r,h_cor,h_overlap]=max_cor(Ch_data,tCor_ST,fsDS,RMG_ch_cor,EMG_ch_cor,OptFig,Chan_Name)



tCor_cal=tCor_ST(1)*fsDS:tCor_ST(2)*fsDS;

tCor_cal=round(tCor_cal');

RMG_temp=Ch_data(tCor_cal,RMG_ch_cor);
EMG_temp=Ch_data(tCor_cal,EMG_ch_cor);
RMG_temp=smoothdata(RMG_temp,'movmean',100);
EMG_temp=smoothdata(EMG_temp,'movmean',100);
RMG_temp=detrend(RMG_temp);
EMG_temp=detrend(EMG_temp);

tOff=((0:(length(RMG_temp)-1))/fsDS)';
[r1,p1] = corrcoef(normalize(RMG_temp),normalize(EMG_temp)); 
r_cor = r1(1,2);
[x,r]=xcorr(RMG_temp,EMG_temp,fsDS,'normalized');  % x+m, y 
max_coef=max(x);
max_r=r(find(x==max_coef))/fsDS;
% if max_r<0 | max_r>0.5
%     RMG_temp=-RMG_temp;
%     fprintf('flip RMG for correlation \n');
%     [r1,p1] = corrcoef(normalize(RMG_temp),normalize(EMG_temp)); 
% r_cor = r1(1,2);
% [x,r]=xcorr(RMG_temp,EMG_temp,fsDS,'normalized');  % x+m, y 
% max_coef=max(x);
% max_r=r(find(x==max_coef))/fsDS;
% end 
%% dtw
% figure()
% dtw(RMG_temp,EMG_temp);


h_cor=[];
h_overlap=[];
if OptFig==1
h_cor=figure();
sz=13;
stem(r/fsDS,x)
hold on
a=linspace(0,1,20);
plot(0*ones(1,length(a)),a,'color',[0.5,0.5,0.5],'LineStyle',':','LineWidth',2);
plot(max_r*ones(1,length(a)),a,'color','r','LineStyle',':','LineWidth',2);
xlabel('Time Lag (s)','FontSize',sz)
ylabel('Cross-Correlation (a.u.)','FontSize',sz)
xlim([-0.5 0.5])
title(['Max Correlation:',num2str(max_coef,'%.2f'),' Time Lag (s):',num2str(max_r,'%.3f')],'FontSize',sz)

set(gcf,'Position',[200,200,600,350]);

h_overlap=figure();
plot(tOff,RMG_temp,'color','red','LineWidth',1.5);
hold on 
plot(tOff,EMG_temp,'color','blue','LineWidth',1.5);
xlabel('Time(s)','FontSize',sz)
ylabel('Amp (a.u.)','FontSize',sz)
title(['Max Correlation:',num2str(max_coef,'%.2f'),' Time Lag (s):',num2str(max_r,'%.3f')],'FontSize',sz)
xlim([0 max(tOff)])
legend({['RMG',' ch',Chan_Name{RMG_ch_cor}],['EMG',' ch',Chan_Name{EMG_ch_cor}]},'FontSize',sz);
set(gcf,'Position',[200,200,2000,300]);
end

end