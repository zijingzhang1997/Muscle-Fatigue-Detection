load(['D:\legRMG\data\','ML_comp_feat\feat_v1_all\feat_v5_all.mat']);
feat_all_final2=feat_all_final;
load(['D:\legRMG\data\','ML_comp_feat\feat_v1_all\feat_v2_all.mat']);
feat_all_final1=feat_all_final;

ind=find(feat_all_final1(20-2,:)~=feat_all_final2(57-2,:));

feat_all_final1(20-3,ind)
feat_all_final2(57-3,ind)

case1=caseName(ind);
