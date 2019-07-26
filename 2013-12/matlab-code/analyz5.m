% // Whole matched set //
clear;

[N,T] = xlsread('forLalit-SupplTable1D-correlations.xlsx');
NSAF = N(:,5:8);
NadjSPC = N(:,9:12);
COV = N(:,14:17);
RPKM = N(:,18:21);

for j=1:4
    fprintf('j=%d: ',j)
    %     if (any(NadjSPC(:,j)==0)), fprintf('NadjSPC '), end
    %     if (any(COV(:,j)==0)), fprintf('COV '), end
    %     if (any(NSAF(:,j)==0)), fprintf('NSAF '), end
    %     if (any(RPKM(:,j)==0)), fprintf('RPKM '), end
    ind = find(NadjSPC(:,j)==0 | COV(:,j)==0);
    if (length(ind)>0), fprintf('NadjSPC-COV %d ',length(ind)), end
    ind = find(NSAF(:,j)==0 | RPKM(:,j)==0);
    if (length(ind)>0), fprintf('NSAF-RPKM %d ',length(ind)), end
    fprintf('\n')
end


fprintf('\n NadjSPC & COV \n');
fprintf(' \t \tPearson\t \tSpearman\n');
fprintf(' \t \tCorrelation\tPvalue\tCorrelation\tPvalue\n');
[Rp,Pp] = corr(log10(mean(NadjSPC,2)),log10(mean(COV,2)),'type','Pearson');
[Rs,Ps] = corr(log10(mean(NadjSPC,2)),log10(mean(COV,2)),'type','Spearman');
fprintf('average\t%d\t%f\t%f\t%f\t%f\n',size(NadjSPC,1),Rp,Pp,Rs,Ps);
sec={'sec1','sec4','sec9','sec14'};
for j=1:4
    ind = find(NadjSPC(:,j)>0 & COV(:,j)>0);
    [Rp,Pp] = corr(log10(NadjSPC(ind,j)),log10(COV(ind,j)),'type','Pearson');
    [Rs,Ps] = corr(log10(NadjSPC(ind,j)),log10(COV(ind,j)),'type','Spearman');
    fprintf('%s\t%d\t%f\t%f\t%f\t%f\n',sec{j},length(ind),Rp,Pp,Rs,Ps);
end

fprintf('\n NSAF & RPKM \n');
fprintf(' \t \tPearson\t \tSpearman\n');
fprintf(' \t \tCorrelation\tPvalue\tCorrelation\tPvalue\n');
[Rp,Pp] = corr(log10(mean(NSAF,2)),log10(mean(RPKM,2)),'type','Pearson');
[Rs,Ps] = corr(log10(mean(NSAF,2)),log10(mean(RPKM,2)),'type','Spearman');
fprintf('average\t%d\t%f\t%f\t%f\t%f\n',size(NadjSPC,1),Rp,Pp,Rs,Ps);
sec={'sec1','sec4','sec9','sec14'};
for j=1:4
    ind = find(NadjSPC(:,j)>0 & COV(:,j)>0);
    [Rp,Pp] = corr(log10(NSAF(ind,j)),log10(RPKM(ind,j)),'type','Pearson');
    [Rs,Ps] = corr(log10(NSAF(ind,j)),log10(RPKM(ind,j)),'type','Spearman');
    fprintf('%s\t%d\t%f\t%f\t%f\t%f\n',sec{j},length(ind),Rp,Pp,Rs,Ps);
end



% fprintf('\n NadjSPC & COV \n');
% fprintf(' \t \tPearson\t \tSpearman\n');
% fprintf(' \t \tCorrelation\tPvalue\tCorrelation\tPvalue\n');
% [Rp,Pp] = corr(log10(mean(NadjSPC,2)),log10(mean(COV,2)),'type','Pearson');
% [Rs,Ps] = corr(log10(mean(NadjSPC,2)),log10(mean(COV,2)),'type','Spearman');
% fprintf('average\t%d\t%f\t%f\t%f\t%f\n',size(NadjSPC,1),Rp,Pp,Rs,Ps);
% sec={'sec1','sec4','sec9','sec14'};
% for j=1:4
%     ind = find(NadjSPC(:,j)>0);
%     [Rp,Pp] = corr(log10(NadjSPC(ind,j)),log10(COV(ind,j)),'type','Pearson');
%     [Rs,Ps] = corr(log10(NadjSPC(ind,j)),log10(COV(ind,j)),'type','Spearman');
%     fprintf('%s\t%d\t%f\t%f\t%f\t%f\n',sec{j},length(ind),Rp,Pp,Rs,Ps);
% end
% 
% fprintf('\n NSAF & RPKM \n');
% fprintf(' \t \tPearson\t \tSpearman\n');
% fprintf(' \t \tCorrelation\tPvalue\tCorrelation\tPvalue\n');
% [Rp,Pp] = corr(log10(mean(NSAF,2)),log10(mean(RPKM,2)),'type','Pearson');
% [Rs,Ps] = corr(log10(mean(NSAF,2)),log10(mean(RPKM,2)),'type','Spearman');
% fprintf('average\t%d\t%f\t%f\t%f\t%f\n',size(NadjSPC,1),Rp,Pp,Rs,Ps);
% sec={'sec1','sec4','sec9','sec14'};
% for j=1:4
%     ind = find(NadjSPC(:,j)>0);
%     [Rp,Pp] = corr(log10(NSAF(ind,j)),log10(RPKM(ind,j)),'type','Pearson');
%     [Rs,Ps] = corr(log10(NSAF(ind,j)),log10(RPKM(ind,j)),'type','Spearman');
%     fprintf('%s\t%d\t%f\t%f\t%f\t%f\n',sec{j},length(ind),Rp,Pp,Rs,Ps);
% end
