function h=plotWaveform_s(Ch_data,ChName,ind,fsDS,titleAll)


tOff=((0:(length(Ch_data)-1))/fsDS)';

h(1)=figure;
a=0;
sz=8;
cN={'red','blue','red','blue','red','blue','red','blue','red','blue','red','blue'};
n=size(Ch_data,2);
for i =1:n
ax(i)=subplot(n,1,i);

plot(tOff,Ch_data(:,i),'LineWidth',0.5,'color',cN{i},'LineWidth',1);

xlabel('time (s)','FontSize',sz)
ylabel('Amp (a.u.)','FontSize',sz)
xlim([a round(max(tOff))])
title([num2str(ind(i)),ChName{i}],'FontSize',sz)


end

set(gcf,'Position',[50,50,1800,150*n]);
sgtitle(titleAll);
linkaxes(ax,'x')

end 