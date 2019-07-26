% // Regression at the gene level //
% already done last time, see sheet5 of 'Supplemental Table 1A-B-C-D-E-aug2013-v1.xls'
% re-doing using log-transformed values
clear;

[N,T] = xlsread('Supplemental Table 1A-B-C-D-E-aug2013-v1.xls',5);
NSAF = N(:,4:7);
RPKM = N(:,8:11);
NadjSPC = N(:,13:16);
COV = N(:,17:20);
groupID = N(:,1);

Y = NSAF; % NadjSPC; % NSAF;
X = RPKM; % COV; % RPKM;

fo=fopen('reg_allnz.txt','w');
fprintf(fo,'GroupID\tIntercept\tSlope\tRsq\tPvalue\n');
% regress only if X is all-nonzero and Y is all-nonzero
for i=1:length(groupID)
    x = X(i,:); y = Y(i,:);
    if any([x, y]==0), continue, end
    stats = regstats(log10(y),log10(x),'linear',{'beta','rsquare','adjrsquare','fstat'});
    fprintf(fo,'%d\t%f\t%f\t%f\t%f\n',groupID(i),stats.beta(1),stats.beta(2),stats.rsquare,stats.fstat.pval);
end
fclose(fo);

fo=fopen('reg_xnz_yrep.txt','w');
fprintf(fo,'GroupID\tIntercept\tSlope\tRsq\tPvalue\n');
% regress only if X is all non-zero and Y zeros replaced with 1e-6
for i=1:length(groupID)
    x = X(i,:); y = Y(i,:);
    if any(x==0), continue, end
    y(find(y==0)) = 1e-6;
    stats = regstats(log10(y),log10(x),'linear',{'beta','rsquare','adjrsquare','fstat'});
    fprintf(fo,'%d\t%f\t%f\t%f\t%f\n',groupID(i),stats.beta(1),stats.beta(2),stats.rsquare,stats.fstat.pval);
end
fclose(fo);

fo=fopen('reg_xrep_yrep.txt','w');
fprintf(fo,'GroupID\tIntercept\tSlope\tRsq\tPvalue\n');
% regress only if X is all non-zero and Y zeros replaced with 1e-6
for i=1:length(groupID)
    x = X(i,:); y = Y(i,:);
    x(find(x==0)) = 1e-6;
    y(find(y==0)) = 1e-6;
    stats = regstats(log10(y),log10(x),'linear',{'beta','rsquare','adjrsquare','fstat'});
    fprintf(fo,'%d\t%f\t%f\t%f\t%f\n',groupID(i),stats.beta(1),stats.beta(2),stats.rsquare,stats.fstat.pval);
end
fclose(fo);

