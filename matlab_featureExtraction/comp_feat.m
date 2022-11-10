%% combine all cases from different experiments 
%% featAllch_all  
%use "genML" script to generate features in cases 
clear all


%fit_parameter_TimeLag=[p(1)*100,p(2),mean(Correct_max_r_seg),mean(delta)];
featAllch_all=[];
fit_parameter_TimeLag_all=[];
label_all=[];
Exp_Name_list={};
File_Name_list={};
NameAll_list={};
fileNameList_0={'Case1fingerInstRoutine1','Case1GripRingRoutine1','Case1gripperRoutine1'...
  };
label_0=[1,1,1];

ExpDate_0='fatigue_Up_v1';

dataPath=['D:\legRMG\data\',ExpDate_0,'\'];
for i=1:length(fileNameList_0)
saveMatFolder=[dataPath,'mat\comp_feat\',fileNameList_0{i}];
load([saveMatFolder,'_Feat_allCh','.mat']);

featAllch_all(:,:,i)=featAll_ch;
fit_parameter_TimeLag_all(:,i)=fit_parameter_TimeLag;
label_all(i)=label_0(i);
Exp_Name_list{i}=[ExpDate_0];
File_Name_list{i}=[fileNameList_0{i}];
NameAll_list{i}=[ExpDate_0,fileNameList_0{i}];
end

num=1;
fileNameList{num}={'Case1fingerInstRoutine2','Case2fingerInstRoutine2',...
    'Case1GripRingRoutine2','Case2GripRingRoutine2',...
    'Case1gripperRoutine2','Case2gripperRoutine2',...
    'Case1bfingerInstRoutine2','Case2bfingerInstRoutine2','Case2bfingerInstRoutine2',...
    'Case1bGripRingRoutine2','Case2bGripRingRoutine2',...
       'Case1p23Routine2','Case2p23Routine2','Case3p23Routine2',...
    'Case1bp23Routine2'};
label_list{num}=[1,1,1,1,1,1,1,1,1,1,1,0,0,0,0];

ExpDate{num}='fatigue_ZZ_hand_new_v1';

%%
% delete 'Case2fatigue_tiptoe_sitRoutine2'
fileNameList{end+1}={'Case1fatigue_tiptoe_standRoutine2','Case1fatigue_tiptoe_sitRoutine2',...
  };
label_list{end+1}=[1,0];

ExpDate{end+1}='fatigue_v1';

%%
% delete 'Case3fatigue_tiptoe_standRoutine1'
fileNameList{end+1}={'Case2fatigue_tiptoe_standRoutine1','Case1fatigue_tiptoe_sitRoutine1','Case2fatigue_tiptoe_sitRoutine1',...
  };
label_list{end+1}=[1,0,0];

ExpDate{end+1}='fatigue_GX';
%% delete 'Case2fatigue_tiptoe_sitRoutine2'
fileNameList{end+1}={'Case1fatigue_tiptoe_standRoutine2'...
  };
label_list{end+1}=[1];

ExpDate{end+1}='fatigue_HS';
%%
fileNameList{end+1}={'Case1fatigue_tiptoe_standRoutine2','Case1fatigue_tiptoe_standRoutine2',...
    'Case1fatigue_tiptoe_sitRoutine2','Case2fatigue_tiptoe_sitRoutine2',...
  };
label_list{end+1}=[1,1,0,0];

ExpDate{end+1}='fatigue_EK';

%%

fileNameList{end+1}={'Case1p23Routine2','Case2p23Routine2',...
    'Case1graspRoutine2',...
    'Case1GripRingRoutine2','Case3GripRingRoutine2',...
   };
label_list{end+1}=[0,0,0,1,1];

ExpDate{end+1}='fatigue_EK_hand_v1';

%%
fileNameList{end+1}={'Case1graspRoutine2',...
    'Case2gripperRoutine2',...
    'Case1fingerInstRoutine2','Case2fingerInstRoutine2',...
    'Case1GripRingRoutine2','Case2GripRingRoutine2',...
   };
label_list{end+1}=[0,1,1,1,1,1];

ExpDate{end+1}='fatigue_KG_hand_v1';

%%
fileNameList{end+1}={'Case1graspRoutine2','Case2graspRoutine2',...
    'Case1p23Routine2',...
    'Case1fingerInstRoutine2','Case2fingerInstRoutine2','Case3fingerInstRoutine2',...
    'Case1GripRingRoutine2','Case2GripRingRoutine2',...
   };
label_list{end+1}=[0,0,0,1,1,1,1,1];

ExpDate{end+1}='fatigue_AK_hand_v1';
%%
fileNameList{end+1}={'Case1graspRoutine2',...
    'Case1p23Routine2','Case2p23Routine2',...
    'Case1fingerInstRoutine2',...
    'Case1GripRingRoutine2',...
   };
label_list{end+1}=[0,0,0,1,1];

ExpDate{end+1}='fatigue_HG_hand_v1';
%%
% TC: 'Case1GripRingRoutine2' =2 so label as non-fatigue
fileNameList{end+1}={'Case1graspRoutine2','Case2graspRoutine2',...
    'Case1p23Routine2','Case2p23Routine2',...
    'Case1fingerInstRoutine2','Case2fingerInstRoutine2',...
    'Case1GripRingRoutine2','Case2GripRingRoutine2',...
    'Case1gripperRoutine2','Case2gripperRoutine2',...
   };
label_list{end+1}=[0,0,0,0,1,1,0,1,1,1];

ExpDate{end+1}='fatigue_TC_hand_v1';
%% import list of experiments 
for comNum=1:length(fileNameList)
    
   
dataPath=['D:\legRMG\data\',ExpDate{comNum},'\'];


for i=1:length(fileNameList{comNum})
saveMatFolder=[dataPath,'mat\comp_feat\',fileNameList{comNum}{i}];
load([saveMatFolder,'_Feat_allCh','.mat']);

featAllch_all(:,:,end+1)=featAll_ch;
fit_parameter_TimeLag_all(:,end+1)=fit_parameter_TimeLag;
label_all(end+1)=label_list{comNum}(i);

Exp_Name_list{end+1}=[ExpDate{comNum}];
File_Name_list{end+1}=[fileNameList{comNum}{i}];
NameAll_list{end+1}=[ExpDate{comNum},fileNameList{comNum}{i}];
end

end
%% all channel feat processing 
% featMean(8) BR PP IN EX mean of (window mean) + BR PP IN EX mean of (window var)
% fit pp: slope *100 (9); p2 (10); delta (11); fit br: slope *100 (12); p2 (13); delta (14);
% err_br : error *100% (15); 
% add *** fit pp: slope/ mean PP (16) *** 

% final feature     15  +4
for i=1: length(label_all)
    featAllch_temp=featAllch_all(:,:,i);
    featAllch_temp(:,end+1) = featAllch_temp(:,9)./featAllch_temp(:,2);
    
    %ind=find(featAllch_temp(:,15)<8);
  if label_all(i)==1
    ind=find(featAllch_temp(:,15)<8 & featAllch_temp(:,9)<0.1);
  end
   if label_all(i)==0
    ind=find(featAllch_temp(:,15)<8);
  end
    %ind=find(featAllch_temp(:,15)<10 & featAllch_temp(:,9)<0);
    if length(ind)==0
        ind=find(featAllch_temp(:,15)<8);
    end
    feat_Passth(i)=length(ind);
feat_th(:,i)=mean(featAllch_temp(ind,:),1);
feat_all_final(:,i)=[feat_th(:,i);fit_parameter_TimeLag_all(:,i)];

end
DataVersionNote='64 cases;  57 featrues; threshold on (15)err_br <8 %; average on all passed channels; label=1 add requirement (9)<0.1 ...fit pp: slope';

save(['D:\legRMG\data\','ML_comp_feat\feat_v1_all\feat_v7_all.mat'],'feat_all_final',...
'label_all','fit_parameter_TimeLag_all','featAllch_all','File_Name_list','Exp_Name_list','DataVersionNote','NameAll_list');
