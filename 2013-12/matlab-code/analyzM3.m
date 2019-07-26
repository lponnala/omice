clear;

FILE='matchedGRMZMdataNoGrouping.txt';
[protID,matched,geneID,cov1,cov4,cov9,cov14,rpkm1,rpkm4,rpkm9,rpkm14,nsaf1,nsaf4,nsaf9,nsaf14,nadjspc1,nadjspc4,nadjspc9,nadjspc14] = textread(FILE,'%s\t%s\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f');
COV=[cov1,cov4,cov9,cov14];
RPKM=[rpkm1,rpkm4,rpkm9,rpkm14];
NSAF=[nsaf1,nsaf4,nsaf9,nsaf14];
NadjSPC=[nadjspc1,nadjspc4,nadjspc9,nadjspc14];

% // correlation using log10 //
% average
fprintf('\n\n--------------------------------\n')
fprintf('average \t %d proteins\n',size(NadjSPC,1));
fprintf('\nlog10(NadjSPC) & log10(COV)\n');
fprintf('Pearson \t \t Spearman \t \n');
fprintf('Rp \t p.value \t Rs \t p.value\n');
[Rp,Pp] = corr(log10(mean(NadjSPC,2)),log10(mean(COV,2)),'type','Pearson');
[Rs,Ps] = corr(log10(mean(NadjSPC,2)),log10(mean(COV,2)),'type','Spearman');
fprintf('%f\t%f\t%f\t%f\n',Rp,Pp,Rs,Ps);
fprintf('\nlog10(NSAF) & log10(RPKM)\n');
fprintf('Pearson \t \t Spearman \t \n');
fprintf('Rp \t p.value \t Rs \t p.value\n');
[Rp,Pp] = corr(log10(mean(NSAF,2)),log10(mean(RPKM,2)),'type','Pearson');
[Rs,Ps] = corr(log10(mean(NSAF,2)),log10(mean(RPKM,2)),'type','Spearman');
fprintf('%f\t%f\t%f\t%f\n',Rp,Pp,Rs,Ps);
% section-wise
sec={'1','4','9','14'};
for j=1:4
    ind = find(NadjSPC(:,j)>0 & COV(:,j)>0 & NSAF(:,j)>0 & RPKM(:,j)>0);
    fprintf('\n\n--------------------------------\n')
    fprintf(['sec',sec{j},' \t %d proteins\n'],length(ind));
    
    fprintf('\nlog10(NadjSPC) & log10(COV)\n');
    fprintf('Pearson \t \t Spearman \t \n');
    fprintf('Rp \t p.value \t Rs \t p.value\n');
    [Rp,Pp] = corr(log10(NadjSPC(ind,j)),log10(COV(ind,j)),'type','Pearson');
    [Rs,Ps] = corr(log10(NadjSPC(ind,j)),log10(COV(ind,j)),'type','Spearman');
    fprintf('%f\t%f\t%f\t%f\n',Rp,Pp,Rs,Ps);
    
    fprintf('\nlog10(NSAF) & log10(RPKM)\n');
    fprintf('Pearson \t \t Spearman \t \n');
    fprintf('Rp \t p.value \t Rs \t p.value\n');
    [Rp,Pp] = corr(log10(NSAF(ind,j)),log10(RPKM(ind,j)),'type','Pearson');
    [Rs,Ps] = corr(log10(NSAF(ind,j)),log10(RPKM(ind,j)),'type','Spearman');
    fprintf('%f\t%f\t%f\t%f\n',Rp,Pp,Rs,Ps);
end
