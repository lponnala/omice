clear;

[N,T] = xlsread('Supplemental Table 1A-B-C-D-E-aug2013-v1.xls',5);
NSAF = N(:,4:7);
RPKM = N(:,8:11);
NadjSPC = N(:,13:16);
COV = N(:,17:20);

% % ---------------------------------------------------
% % // PLOT DISTRIBUTIONS //
% % ---------------------------------------------------
% colors = ['r','g','b','m'];
% figure;
% for j=1:4
%     [f,xi]=ksdensity(NadjSPC(:,j));
%     plot(xi,f,colors(j))
%     hold on
% end
% xlabel('NadjSPC')
% legend('sec1','sec4','sec9','sec14','Location','NorthEast')
% hold off
% saveas(gcf,'NadjSPC','jpg')
% figure;
% for j=1:4
%     [f,xi]=ksdensity(COV(:,j));
%     plot(xi,f,colors(j))
%     hold on
% end
% xlabel('COV')
% legend('sec1','sec4','sec9','sec14','Location','NorthEast')
% hold off
% saveas(gcf,'COV','jpg')
% figure;
% for j=1:4
%     [f,xi]=ksdensity(NSAF(:,j));
%     plot(xi,f,colors(j))
%     hold on
% end
% xlabel('NSAF')
% legend('sec1','sec4','sec9','sec14','Location','NorthEast')
% hold off
% saveas(gcf,'NSAF','jpg')
% figure;
% for j=1:4
%     [f,xi]=ksdensity(RPKM(:,j));
%     plot(xi,f,colors(j))
%     hold on
% end
% xlabel('RPKM')
% legend('sec1','sec4','sec9','sec14','Location','NorthEast')
% hold off
% saveas(gcf,'RPKM','jpg')
% % ---------------------------------------------------


% ---------------------------------------------------
% // CALCULATE CORRELATIONS & MAKE PLOTS //
% ---------------------------------------------------

cutoff = 0; % 0.0002; 0.00005;
sec = {'sec1','sec4','sec9','sec14'};

fprintf('\n NadjSPC & COV \n')
% average
NadjSPC_avg = mean(NadjSPC,2); COV_avg = mean(COV,2);
ind = find(NadjSPC_avg>=cutoff);
plot(log10(NadjSPC_avg(ind)),log10(COV_avg(ind)),'b*');
title(['Using proteins whose avg NadjSPC >= ',num2str(cutoff)]); xlabel('log10(NadjSPC)'); ylabel('log10(COV)');
saveas(gcf, ['avg_NadjSPC_COV_',num2str(cutoff),'.jpg'], 'jpg');
[Rp,Pp] = corr(log10(NadjSPC_avg(ind)),log10(COV_avg(ind)),'type','Pearson');
[Rs,Ps] = corr(log10(NadjSPC_avg(ind)),log10(COV_avg(ind)),'type','Spearman');
fprintf('average\t%d\t%f\t%f\t%f\t%f\n',length(ind),Rp,Pp,Rs,Ps);
% section-wise
ind = cell(1,4);
for j=1:4
    ind{j} = find(NadjSPC(:,j)>=cutoff);
    [Rp,Pp] = corr(log10(NadjSPC(ind{j},j)),log10(COV(ind{j},j)),'type','Pearson');
    [Rs,Ps] = corr(log10(NadjSPC(ind{j},j)),log10(COV(ind{j},j)),'type','Spearman');
    fprintf('%s\t%d\t%f\t%f\t%f\t%f\n',sec{j},length(ind{j}),Rp,Pp,Rs,Ps);
end
plot(log10(NadjSPC(ind{1},1)),log10(COV(ind{1},1)),'r.',log10(NadjSPC(ind{2},2)),log10(COV(ind{2},2)),'g.',log10(NadjSPC(ind{3},3)),log10(COV(ind{3},3)),'b.',log10(NadjSPC(ind{4},4)),log10(COV(ind{4},4)),'m.')
title(['Using proteins for which NadjSPC >= ',num2str(cutoff)]); xlabel('log10(NadjSPC)'); ylabel('log10(COV)');
legend('sec1','sec4','sec9','sec14','Location','NorthWest')
saveas(gcf,['sec_NadjSPC_',num2str(cutoff),'.jpg'],'jpg')

fprintf('\n NSAF & RPKM \n')
% average
NSAF_avg = mean(NSAF,2); RPKM_avg = mean(RPKM,2);
ind = find(NadjSPC_avg>=cutoff);
plot(log10(NSAF_avg(ind)),log10(RPKM_avg(ind)),'b*');
title(['Using proteins whose avg NSAF >= ',num2str(cutoff)]); xlabel('log10(NSAF)'); ylabel('log10(RPKM)');
saveas(gcf, ['avg_NSAF_RPKM_',num2str(cutoff),'.jpg'], 'jpg');
[Rp,Pp] = corr(log10(NSAF_avg(ind)),log10(RPKM_avg(ind)),'type','Pearson');
[Rs,Ps] = corr(log10(NSAF_avg(ind)),log10(RPKM_avg(ind)),'type','Spearman');
fprintf('average\t%d\t%f\t%f\t%f\t%f\n',length(ind),Rp,Pp,Rs,Ps);
% section-wise
ind = cell(1,4);
for j=1:4
    ind{j} = find(NadjSPC(:,j)>=cutoff);
    [Rp,Pp] = corr(log10(NSAF(ind{j},j)),log10(RPKM(ind{j},j)),'type','Pearson');
    [Rs,Ps] = corr(log10(NSAF(ind{j},j)),log10(RPKM(ind{j},j)),'type','Spearman');
    fprintf('%s\t%d\t%f\t%f\t%f\t%f\n',sec{j},length(ind{j}),Rp,Pp,Rs,Ps);
end
plot(log10(NSAF(ind{1},1)),log10(RPKM(ind{1},1)),'r.',log10(NSAF(ind{2},2)),log10(RPKM(ind{2},2)),'g.',log10(NSAF(ind{3},3)),log10(RPKM(ind{3},3)),'b.',log10(NSAF(ind{4},4)),log10(RPKM(ind{4},4)),'m.')
title(['Using proteins for which NSAF >= ',num2str(cutoff)]); xlabel('log10(NSAF)'); ylabel('log10(RPKM)');
legend('sec1','sec4','sec9','sec14','Location','NorthWest')
saveas(gcf,['sec_NSAF_',num2str(cutoff),'.jpg'],'jpg')
% ---------------------------------------------------
