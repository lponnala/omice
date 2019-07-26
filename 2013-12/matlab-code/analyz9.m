% // Clustering protein/mRNA ratio at the gene level //
clear;

[N,T] = xlsread('Supplemental Table 1A-B-C-D-E-aug2013-v1.xls',5);
NSAF = N(:,4:7);
RPKM = N(:,8:11);
NadjSPC = N(:,13:16);
COV = N(:,17:20);
groupID = N(:,1);

num = NSAF; % NadjSPC;
den = RPKM; % COV;

% don't replace 0 values of num
repsmall = 0;

% % replace 0 values of num with 1e-6
% repsmall = 1;
% num(find(num==0))=1e-6;

ratio = num./den;
cols = {'sec1','sec4','sec9','sec14'};

% group genes by assigned cluster
if repsmall
    [N,T]=xlsread('NSAF_RPKM_repsmall.xlsx');
else
    [N,T]=xlsread('NSAF_RPKM.xlsx');
end
numa = nan(size(N,2),length(cols)); numm = nan(size(N,2),length(cols));
dena = nan(size(N,2),length(cols)); denm = nan(size(N,2),length(cols));
ratioa = nan(size(N,2),length(cols)); ratiom = nan(size(N,2),length(cols));
for j=1:size(N,2)
    ids = N(2:end,j);
    ids = ids(find(isfinite(ids)));
    cr=[]; cn=[]; cd=[];
    for i=1:length(ids)
        ind = find(groupID==ids(i));
        if (length(ind)~=1), error('ind not found!'), end
        cr = [cr; ratio(ind,:)];
        cn = [cn; num(ind,:)];
        cd = [cd; den(ind,:)];
    end
    figure
    hold on, plot(1:4,mean(cr),'r*'), plot(1:4,median(cr),'b*'), hold off
    axis([0 5 min([mean(cr) median(cr)]) 1.1*max([mean(cr) median(cr)])])
    set(gca,'XTick',[1 2 3 4],'XTickLabel',cols)
    title(['cluster',num2str(j),' : NSAF/RPKM ratio'])
    legend('average','median','Location','Best')
    if repsmall
        saveas(gcf,['nsaf_rpkm_repsmall.cluster_',num2str(j),'.jpg'],'jpg')
    else
        saveas(gcf,['nsaf_rpkm.cluster_',num2str(j),'.jpg'],'jpg')
    end
    
    numa(j,:) = mean(cn); numm(j,:) = median(cn);
    dena(j,:) = mean(cd); denm(j,:) = median(cd);
    ratioa(j,:) = mean(cr); ratiom(j,:) = median(cr);
end

% print average & median values of ratio to a spreadsheet
if repsmall
    xlswrite('nsaf_rpkm_repsmall.ratio_stats.xls',ratioa,1);
    xlswrite('nsaf_rpkm_repsmall.ratio_stats.xls',ratiom,2);
else
    xlswrite('nsaf_rpkm.ratio_stats.xls',ratioa,1);
    xlswrite('nsaf_rpkm.ratio_stats.xls',ratiom,2);
end

symb = {'r*','b*','g*','c*','m*','k*'};

% plot average values of mrna-prot abundance
leg = {};
figure;
hold on
for i=1:size(numa,1)
    plot(log10(dena(i,:)),log10(numa(i,:)),symb{i})
    leg = [leg, ['cluster',num2str(i)]];
end
hold off
xlabel('log10(RPKM)'); ylabel('log10(NSAF)');
legend(leg,'Location','Best')
if repsmall
    saveas(gcf,'nsaf_rpkm_repsmall.cluster_average.jpg','jpg')
else
    saveas(gcf,'nsaf_rpkm.cluster_average.jpg','jpg')
end

% plot median values of mrna-prot abundance
leg = {};
figure;
hold on
for i=1:size(numm,1)
    plot(log10(denm(i,:)),log10(numm(i,:)),symb{i})
    leg = [leg, ['cluster',num2str(i)]];
end
hold off
xlabel('log10(RPKM)'); ylabel('log10(NSAF)');
legend(leg,'Location','Best')
if repsmall
    saveas(gcf,'nsaf_rpkm_repsmall.cluster_median.jpg','jpg')
else
    saveas(gcf,'nsaf_rpkm.cluster_median.jpg','jpg')
end
