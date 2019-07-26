% // 70S & 80S //
clear;

% ----------------------------------------------------------------------
% "function x section" average abundances
% ----------------------------------------------------------------------
[N,T] = xlsread('forLalit-70Svs80S-correlations.xlsx');
NSAF = N(:,11:14);
RPKM = N(:,15:18);
NadjSPC = N(:,20:23);
COV = N(:,24:27);

func = unique(T(2:end,7));
sec = {'sec1','sec4','sec9','sec14'};

fprintf('\n Using all proteins \n');
fprintf('Section\t#Proteins\tavg_NadjSPC\tavg_COV\tavg_NSAF\tavg_RPKM\n');
for i=1:length(func)
    fprintf('Function: %s\n',func{i});
    for j=1:length(sec)
        ind = find(strcmp(T(2:end,7),func{i}));
        fprintf('%s\t%d\t%f\t%f\t%f\t%f\n',sec{j},length(ind),mean(NadjSPC(ind,j)),mean(COV(ind,j)),mean(NSAF(ind,j)),mean(RPKM(ind,j)));
    end
end

fprintf('\n Excluding proteins that have NadjSPC=0 \n');
fprintf('Section\t#Proteins\tavg_NadjSPC\tavg_COV\tavg_NSAF\tavg_RPKM\n');
for i=1:length(func)
    fprintf('Function: %s\n',func{i});
    for j=1:length(sec)
        ind = find(strcmp(T(2:end,7),func{i}) & NadjSPC(:,j)>0);
        fprintf('%s\t%d\t%f\t%f\t%f\t%f\n',sec{j},length(ind),mean(NadjSPC(ind,j)),mean(COV(ind,j)),mean(NSAF(ind,j)),mean(RPKM(ind,j)));
    end
end
% ----------------------------------------------------------------------


% % ----------------------------------------------------------------------
% % Within each section, leave out zero values, with log transformation
% % ----------------------------------------------------------------------
% [N,T] = xlsread('forLalit-70Svs80S-correlations.xlsx');
% NSAF = N(:,11:14);
% RPKM = N(:,15:18);
% NadjSPC = N(:,20:23);
% COV = N(:,24:27);
% 
% func = unique(T(2:end,7));
% sec = {'sec1','sec4','sec9','sec14'};
% 
% fprintf('\n NadjSPC & COV \n');
% fprintf('Section\t#Proteins\tRp\tp.val\tRs\tp.val\n');
% for i=1:length(func)
%     fprintf('Function: %s\n',func{i});
%     for j=1:length(sec)
%         ind = find(strcmp(T(2:end,7),func{i}) & NadjSPC(:,j)>0);
%         [Rp,Pp] = corr(log10(NadjSPC(ind,j)),log10(COV(ind,j)),'type','Pearson');
%         [Rs,Ps] = corr(log10(NadjSPC(ind,j)),log10(COV(ind,j)),'type','Spearman');
%         fprintf('%s\t%d\t%f\t%f\t%f\t%f\n',sec{j},length(ind),Rp,Pp,Rs,Ps);
%     end
% end
% 
% fprintf('\n NSAF & RPKM \n');
% fprintf('Section\t#Proteins\tRp\tp.val\tRs\tp.val\n');
% for i=1:length(func)
%     fprintf('Function: %s\n',func{i});
%     for j=1:length(sec)
%         ind = find(strcmp(T(2:end,7),func{i}) & NadjSPC(:,j)>0);
%         [Rp,Pp] = corr(log10(NSAF(ind,j)),log10(RPKM(ind,j)),'type','Pearson');
%         [Rs,Ps] = corr(log10(NSAF(ind,j)),log10(RPKM(ind,j)),'type','Spearman');
%         fprintf('%s\t%d\t%f\t%f\t%f\t%f\n',sec{j},length(ind),Rp,Pp,Rs,Ps);
%     end
% end
% % ----------------------------------------------------------------------


% % ----------------------------------------------------------------------
% % Within each section, leave out zero values, no log transformation
% % ----------------------------------------------------------------------
% [N,T] = xlsread('forLalit-70Svs80S-correlations.xlsx');
% NSAF = N(:,11:14);
% RPKM = N(:,15:18);
% NadjSPC = N(:,20:23);
% COV = N(:,24:27);
% 
% func = unique(T(2:end,7));
% sec = {'sec1','sec4','sec9','sec14'};
% 
% fprintf('\n NadjSPC & COV \n');
% fprintf('Section\t#Proteins\tRp\tp.val\tRs\tp.val\n');
% for i=1:length(func)
%     fprintf('Function: %s\n',func{i});
%     for j=1:length(sec)
%         ind = find(strcmp(T(2:end,7),func{i}) & NadjSPC(:,j)>0);
%         [Rp,Pp] = corr(NadjSPC(ind,j),COV(ind,j),'type','Pearson');
%         [Rs,Ps] = corr(NadjSPC(ind,j),COV(ind,j),'type','Spearman');
%         fprintf('%s\t%d\t%f\t%f\t%f\t%f\n',sec{j},length(ind),Rp,Pp,Rs,Ps);
%     end
% end
% 
% fprintf('\n NSAF & RPKM \n');
% fprintf('Section\t#Proteins\tRp\tp.val\tRs\tp.val\n');
% for i=1:length(func)
%     fprintf('Function: %s\n',func{i});
%     for j=1:length(sec)
%         ind = find(strcmp(T(2:end,7),func{i}) & NadjSPC(:,j)>0);
%         [Rp,Pp] = corr(NSAF(ind,j),RPKM(ind,j),'type','Pearson');
%         [Rs,Ps] = corr(NSAF(ind,j),RPKM(ind,j),'type','Spearman');
%         fprintf('%s\t%d\t%f\t%f\t%f\t%f\n',sec{j},length(ind),Rp,Pp,Rs,Ps);
%     end
% end
% % ----------------------------------------------------------------------


% % ----------------------------------------------------------------------
% % Within each section, no cutoff
% % ----------------------------------------------------------------------
% [N,T] = xlsread('forLalit-70Svs80S-correlations.xlsx');
% NSAF = N(:,11:14);
% RPKM = N(:,15:18);
% NadjSPC = N(:,20:23);
% COV = N(:,24:27);
% 
% func = unique(T(2:end,7));
% sec = {'sec1','sec4','sec9','sec14'};
% 
% fprintf('\n NadjSPC & COV \n');
% fprintf('Section\t#Proteins\tRp\tp.val\tRs\tp.val\n');
% for i=1:length(func)
%     ind = find(strcmp(T(2:end,7),func{i}));
%     fprintf('Function: %s\n',func{i});
%     for j=1:length(sec)
%         [Rp,Pp] = corr(NadjSPC(ind,j),COV(ind,j),'type','Pearson');
%         [Rs,Ps] = corr(NadjSPC(ind,j),COV(ind,j),'type','Spearman');
%         fprintf('%s\t%d\t%f\t%f\t%f\t%f\n',sec{j},length(ind),Rp,Pp,Rs,Ps);
%     end
% end
% 
% fprintf('\n NSAF & RPKM \n');
% fprintf('Section\t#Proteins\tRp\tp.val\tRs\tp.val\n');
% for i=1:length(func)
%     ind = find(strcmp(T(2:end,7),func{i}));
%     fprintf('Function: %s\n',func{i});
%     for j=1:length(sec)
%         [Rp,Pp] = corr(NSAF(ind,j),RPKM(ind,j),'type','Pearson');
%         [Rs,Ps] = corr(NSAF(ind,j),RPKM(ind,j),'type','Spearman');
%         fprintf('%s\t%d\t%f\t%f\t%f\t%f\n',sec{j},length(ind),Rp,Pp,Rs,Ps);
%     end
% end
% % ----------------------------------------------------------------------


% % ----------------------------------------------------------------------
% % Using average across sections
% % ----------------------------------------------------------------------
% [N,T] = xlsread('forLalit-70Svs80S-correlations.xlsx');
% NSAF = N(:,11:14);
% RPKM = N(:,15:18);
% NadjSPC = N(:,20:23);
% COV = N(:,24:27);
% 
% func = unique(T(2:end,7));
% 
% fprintf('\n NadjSPC & COV \n');
% fprintf(' \t \tPearson\t \tSpearman\n');
% fprintf(' \t \tCorrelation\tPvalue\tCorrelation\tPvalue\n');
% NadjSPC_avg = mean(NadjSPC,2); COV_avg = mean(COV,2);
% for j=1:length(func)
%     ind = find(strcmp(T(2:end,7),func{j}));
%     [Rp,Pp] = corr(log10(NadjSPC_avg(ind)),log10(COV_avg(ind)),'type','Pearson');
%     [Rs,Ps] = corr(log10(NadjSPC_avg(ind)),log10(COV_avg(ind)),'type','Spearman');
%     fprintf('%s\t%d\t%f\t%f\t%f\t%f\n',func{j},length(ind),Rp,Pp,Rs,Ps);
% end
% 
% fprintf('\n NSAF & RPKM \n');
% fprintf(' \t \tPearson\t \tSpearman\n');
% fprintf(' \t \tCorrelation\tPvalue\tCorrelation\tPvalue\n');
% NSAF_avg = mean(NSAF,2); RPKM_avg = mean(RPKM,2);
% for j=1:length(func)
%     ind = find(strcmp(T(2:end,7),func{j}));
%     [Rp,Pp] = corr(log10(NSAF_avg(ind)),log10(RPKM_avg(ind)),'type','Pearson');
%     [Rs,Ps] = corr(log10(NSAF_avg(ind)),log10(RPKM_avg(ind)),'type','Spearman');
%     fprintf('%s\t%d\t%f\t%f\t%f\t%f\n',func{j},length(ind),Rp,Pp,Rs,Ps);
% end
% % ----------------------------------------------------------------------
