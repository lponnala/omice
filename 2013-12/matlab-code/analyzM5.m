clear;

[N,T] = xlsread('matchedGRMZMdataNoGrouping.xlsx');
COV=N(:,1:4);
RPKM=N(:,5:8);
groupID=N(:,9);
NSAF=N(:,10:13);
NadjSPC=N(:,14:17);

% look for number of group IDs among non-zero proteins in each section
sec = {'sec1','sec4','sec9','sec14'};

fprintf('\nGROUPED PROTEINS\n');
fprintf('section\tProteins\tGroups\tCOV\t\t\tRPKM\t\t\tNSAF\t\t\tNadjSPC\t\t\n');
fprintf('\t\t\tmean\tmedian\tstdev\tmean\tmedian\tstdev\tmean\tmedian\tstdev\tmean\tmedian\tstdev\n');
for j=1:4
    % gi = find(NadjSPC(:,j)>0 & COV(:,j)>0 & NSAF(:,j)>0 & RPKM(:,j)>0 & isfinite(groupID));
    % gi = find(isfinite(groupID));
    gi = find(NadjSPC(:,j)>0 & isfinite(groupID));
    
    fprintf('%s\t%d\t%d\t',sec{j},length(gi),length(unique(groupID(gi))));
    fprintf('%f\t%f\t%f\t',mean(COV(gi,j)),median(COV(gi,j)),std(COV(gi,j)));
    fprintf('%f\t%f\t%f\t',mean(RPKM(gi,j)),median(RPKM(gi,j)),std(RPKM(gi,j)));
    fprintf('%f\t%f\t%f\t',mean(NSAF(gi,j)),median(NSAF(gi,j)),std(NSAF(gi,j)));
    fprintf('%f\t%f\t%f\n',mean(NadjSPC(gi,j)),median(NadjSPC(gi,j)),std(NadjSPC(gi,j)));
end

fprintf('\nUNGROUPED PROTEINS\n');
fprintf('section\tProteins\tGroups\tCOV\t\t\tRPKM\t\t\tNSAF\t\t\tNadjSPC\t\t\n');
fprintf('\t\t\tmean\tmedian\tstdev\tmean\tmedian\tstdev\tmean\tmedian\tstdev\tmean\tmedian\tstdev\n');
for j=1:4
    % ui = find(NadjSPC(:,j)>0 & COV(:,j)>0 & NSAF(:,j)>0 & RPKM(:,j)>0 & ~isfinite(groupID));
    % ui = find(~isfinite(groupID));
    ui = find(NadjSPC(:,j)>0 & ~isfinite(groupID));
    
    fprintf('%s\t%d\t%d\t',sec{j},length(ui),length(unique(groupID(ui))));
    fprintf('%f\t%f\t%f\t',mean(COV(ui,j)),median(COV(ui,j)),std(COV(ui,j)));
    fprintf('%f\t%f\t%f\t',mean(RPKM(ui,j)),median(RPKM(ui,j)),std(RPKM(ui,j)));
    fprintf('%f\t%f\t%f\t',mean(NSAF(ui,j)),median(NSAF(ui,j)),std(NSAF(ui,j)));
    fprintf('%f\t%f\t%f\n',mean(NadjSPC(ui,j)),median(NadjSPC(ui,j)),std(NadjSPC(ui,j)));
end

% plot the data
[N,T]=xlsread('gStat_AllData.xlsx');
sec = 1:4;
R = [];
r = N(6:9,1)./N(1:4,1); R = [R; r'];
r = N(6:9,4)./N(1:4,4); R = [R; r'];
r = N(6:9,7)./N(1:4,7); R = [R; r'];
r = N(6:9,10)./N(1:4,10); R = [R; r'];
plot(sec,R(1,:),'r*',sec,R(2,:),'g*')
legend('COV','RPKM')
set(gca,'XTick',sec)
set(gca,'XTickLabel',{'sec1','sec4','sec9','sec14'})


