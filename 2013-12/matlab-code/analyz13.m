clear;

[N,T] = xlsread('calculate-correlation.xlsx');

% NSAF = N(1:2:22,3:6);
% RPKM = N(2:2:22,3:6);

NSAF = N(1:2:18,3:6);
RPKM = N(2:2:18,3:6);

fprintf('\n -- using raw values --\n');
[Rp,Pp] = corr(NSAF(:),RPKM(:),'type','Pearson');
[Rs,Ps] = corr(NSAF(:),RPKM(:),'type','Spearman');
fprintf('Type\tCorrelation\tP-value\n');
fprintf('Pearson (Rp)\t%f\t%f\n',Rp,Pp);
fprintf('Spearman (Rs)\t%f\t%f\n',Rs,Ps);

fprintf('\n -- using log10 values --\n');
[Rp,Pp] = corr(log10(NSAF(:)),log10(RPKM(:)),'type','Pearson');
[Rs,Ps] = corr(log10(NSAF(:)),log10(RPKM(:)),'type','Spearman');
fprintf('Type\tCorrelation\tP-value\n');
fprintf('Pearson (Rp)\t%f\t%f\n',Rp,Pp);
fprintf('Spearman (Rs)\t%f\t%f\n',Rs,Ps);
