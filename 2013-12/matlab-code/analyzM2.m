clear;

FILE='matchedGRMZMdataNoGrouping.txt';
[protID,matched,geneID,cov1,cov4,cov9,cov14,rpkm1,rpkm4,rpkm9,rpkm14,nsaf1,nsaf4,nsaf9,nsaf14,nadjspc1,nadjspc4,nadjspc9,nadjspc14] = textread(FILE,'%s\t%s\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f');
COV=[cov1,cov4,cov9,cov14];
RPKM=[rpkm1,rpkm4,rpkm9,rpkm14];
NSAF=[nsaf1,nsaf4,nsaf9,nsaf14];
NadjSPC=[nadjspc1,nadjspc4,nadjspc9,nadjspc14];

% // correlation using log10 //
% average
[Rp,Pp] = corr(log10(mean(NSAF,2)),log10(mean(RPKM,2)),'type','Pearson');
[Rs,Ps] = corr(log10(mean(NSAF,2)),log10(mean(RPKM,2)),'type','Spearman');
fprintf('NSAF & RPKM\n');
fprintf('Pearson (Rp)\t%f\t%f\n',Rp,Pp);
fprintf('Spearman (Rs)\t%f\t%f\n',Rs,Ps);
[Rp,Pp] = corr(log10(mean(NadjSPC,2)),log10(mean(COV,2)),'type','Pearson');
[Rs,Ps] = corr(log10(mean(NadjSPC,2)),log10(mean(COV,2)),'type','Spearman');
fprintf('NadjSPC & COV\n');
fprintf('Pearson (Rp)\t%f\t%f\n',Rp,Pp);
fprintf('Spearman (Rs)\t%f\t%f\n',Rs,Ps);
% section-wise
sec={'1','4','9','14'};
COV(find(COV==0))=1e-6;
RPKM(find(RPKM==0))=1e-6;
NSAF(find(NSAF==0))=1e-6;
NadjSPC(find(NadjSPC==0))=1e-6;
for j=1:4
    [Rp,Pp] = corr(log10(NSAF(:,j)),log10(RPKM(:,j)),'type','Pearson');
    [Rs,Ps] = corr(log10(NSAF(:,j)),log10(RPKM(:,j)),'type','Spearman');
    fprintf(['sec',sec{j},' : NSAF & RPKM\n']);
    fprintf('Pearson (Rp)\t%f\t%f\n',Rp,Pp);
    fprintf('Spearman (Rs)\t%f\t%f\n',Rs,Ps);
    [Rp,Pp] = corr(log10(NadjSPC(:,j)),log10(COV(:,j)),'type','Pearson');
    [Rs,Ps] = corr(log10(NadjSPC(:,j)),log10(COV(:,j)),'type','Spearman');
    fprintf(['sec',sec{j},' : NadjSPC & COV\n']);
    fprintf('Pearson (Rp)\t%f\t%f\n',Rp,Pp);
    fprintf('Spearman (Rs)\t%f\t%f\n',Rs,Ps);
end
