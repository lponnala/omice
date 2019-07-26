clear;

[N,T] = xlsread('matchedGRMZMdataNoGrouping.xlsx');
COV=N(:,1:4);
RPKM=N(:,5:8);
groupID=N(:,9);
NSAF=N(:,10:13);
NadjSPC=N(:,14:17);

% look for number of group IDs among non-zero proteins in each section
sec = {'sec1','sec4','sec9','sec14'};
fprintf('section\tProteins (P)\tGroupedProteins (GP)\tGroups (G)\n');
for j=1:4
    ind = find(NadjSPC(:,j)>0 & COV(:,j)>0 & NSAF(:,j)>0 & RPKM(:,j)>0);
    gid = groupID(ind);
    fprintf('%s\t%d\t%d\t%d\n',sec{j},length(ind),sum(isfinite(groupID(ind))),length(unique(gid(isfinite(gid)))));
end

