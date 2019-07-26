clear;

cols={'W0-1','W2-3','W4-5','W8-9','W12-13','M0-1','M2-3','M4-5','M8-9','M12-13'};
outs = {'WT_MUT_20','WT_MUT_10'};
for k=1:length(outs)
	file=[outs{k},'.txt'];
	[ID,w1,w2,w3,w4,w5,m1,m2,m3,m4,m5]=textread(file,'%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f');
	Data=[w1 w2 w3 w4 w5 m1 m2 m3 m4 m5];
	outname=outs{k};
	CGobj=clustergram(Data,'RowLabels',ID,'ColumnLabels',cols,'RowPdist','correlation','Linkage','average','Cluster',1);
	plot(CGobj); title([outname,' clusters']); saveas(gcf,[outname,'.jpg'],'jpg');
	save([outname,'_obj.mat'],'CGobj');
	clear file ID w1 w2 w3 w4 w5 m1 m2 m3 m4 m5 CGobj
end


cols={'0-1','2-3','4-5','8-9','12-13'};
outs = {'WT_20','WT_10'};
for k=1:length(outs)
	file=[outs{k},'.txt'];
	[ID,w1,w2,w3,w4,w5]=textread(file,'%s\t%f\t%f\t%f\t%f\t%f');
	Data=[w1 w2 w3 w4 w5];
	outname=outs{k};
	CGobj=clustergram(Data,'RowLabels',ID,'ColumnLabels',cols,'RowPdist','correlation','Linkage','average','Cluster',1);
	plot(CGobj); title([outname,' clusters']); saveas(gcf,[outname,'.jpg'],'jpg');
	save([outname,'_obj.mat'],'CGobj');
	clear file ID w1 w2 w3 w4 w5 CGobj
end


% [ID,w1,w2,w3,w4,w5]=textread('WT_10.txt','%s\t%f\t%f\t%f\t%f\t%f');
% Data=[w1 w2 w3 w4 w5];
% outname='WT_10';
% CGobj=clustergram(Data,'RowLabels',ID,'ColumnLabels',cols,'RowPdist','correlation','Linkage','average','Cluster',1);
% plot(CGobj); title([outname,' clusters']); saveas(gcf,[outname,'.jpg'],'jpg');
% save([outname,'_obj.mat'],'CGobj');


% [ID,w1,w2,w3,w4,w5,m1,m2,m3,m4,m5]=textread('g10.txt','%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f');
% cols={'0-1','2-3','4-5','8-9','12-13'};
% 
% W=[w1 w2 w3 w4 w5];
% outname='WT_10';
% CGobj=clustergram(W,'RowLabels',ID,'ColumnLabels',cols,'RowPdist','correlation','Linkage','average','Cluster',1);
% % CGobj=clustergram(Data,'RowLabels',ID,'ColumnLabels',cols,'RowPdist','spearman','Linkage','average','Cluster',1);
% plot(CGobj); title([outname,' clusters']); saveas(gcf,[outname,'.jpg'],'jpg');
% save([outname,'_obj.mat'],'CGobj');
% 
% M=[m1 m2 m3 m4 m5];
% outname='mutant_10';
% CGobj=clustergram(M,'RowLabels',ID,'ColumnLabels',cols,'RowPdist','correlation','Linkage','average','Cluster',1);
% % CGobj=clustergram(Data,'RowLabels',ID,'ColumnLabels',cols,'RowPdist','spearman','Linkage','average','Cluster',1);
% plot(CGobj); title([outname,' clusters']); saveas(gcf,[outname,'.jpg'],'jpg');
% save([outname,'_obj.mat'],'CGobj');
