clear;

FILE = 'allData.txt';
[protID,matched,geneID,cov1,cov4,cov9,cov14,rpkm1,rpkm4,rpkm9,rpkm14,groupID,nsaf1,nsaf4,nsaf9,nsaf14,nadjspc1,nadjspc4,nadjspc9,nadjspc14] = textread(FILE,'%s\t%s\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f');
COV=[cov1,cov4,cov9,cov14];
RPKM=[rpkm1,rpkm4,rpkm9,rpkm14];
NSAF=[nsaf1,nsaf4,nsaf9,nsaf14];
NadjSPC=[nadjspc1,nadjspc4,nadjspc9,nadjspc14];

% // cross correlation plots //
% average
figure, plot(log10(mean(NSAF,2)),log10(mean(RPKM,2)),'b*'); xlabel('log10(NSAF)'); ylabel('log10(RPKM)'); title('Using 113 matched genes, average across all sections'); saveas(gcf,['nsaf_rpkm_avg.jpg'],'jpg');
figure, plot(log10(mean(NadjSPC,2)),log10(mean(COV,2)),'b*'); xlabel('log10(NadjSPC)'); ylabel('log10(COV)'); title('Using 113 matched genes, average across all sections'); saveas(gcf,['nadjspc_cov_avg.jpg'],'jpg');
% % section-wise
% sec={'1','4','9','14'};
% for j=1:4
%     figure, plot(log10(NSAF(:,j)),log10(RPKM(:,j)),'b*'); xlabel('log10(NSAF)'); ylabel('log10(RPKM)'); title(['Using 113 matched genes, section ',sec{j}]); saveas(gcf,['nsaf_rpkm_',sec{j},'.jpg'],'jpg');
%     figure, plot(log10(NadjSPC(:,j)),log10(COV(:,j)),'b*'); xlabel('log10(NadjSPC)'); ylabel('log10(COV)'); title(['Using 113 matched genes, section ',sec{j}]); saveas(gcf,['nadjspc_cov_',sec{j},'.jpg'],'jpg');
% end

% // correlation //
% average, log10
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
