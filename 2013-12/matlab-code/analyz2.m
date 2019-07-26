% // Plastid or non-plastid //
clear;

[N,T] = xlsread('forLalit-SupplTable1E-by-plastid.xlsx');
NSAF = N(:,8:11);
RPKM = N(:,12:15);
NadjSPC = N(:,17:20);
COV = N(:,21:24);
ptd = T(2:end,3);

cutoff = 0.00005;
sec = {'sec1','sec4','sec9','sec14'};

fprintf('\n NadjSPC & COV \n')
% -- plastid genes --
fprintf('\n-- Plastid proteins --\n');
fprintf('\nCORRELATION USING LOG10\n')
fprintf('type\t#prot\tRp\tp.val\tRs\tp.val\n')
ind_p = find(strcmp(ptd,'yes'));
[Rp,Pp] = corr(log10(mean(NadjSPC(ind_p,:),2)), log10(mean(COV(ind_p,:),2)), 'type', 'Pearson');
[Rs,Ps] = corr(log10(mean(NadjSPC(ind_p,:),2)), log10(mean(COV(ind_p,:),2)), 'type', 'Spearman');
fprintf('average\t%d\t%f\t%f\t%f\t%f\n',length(ind_p),Rp,Pp,Rs,Ps);
for j=1:4
    ind = find(strcmp(ptd,'yes') & NadjSPC(:,j)>=cutoff);
    [rp,pp] = corr(log10(NadjSPC(ind,j)),log10(COV(ind,j)),'type','Pearson');
    [rs,ps] = corr(log10(NadjSPC(ind,j)),log10(COV(ind,j)),'type','Spearman');
    fprintf('%s\t%d\t%f\t%f\t%f\t%f\n',sec{j},length(ind),rp,pp,rs,ps);
end
fprintf('\nAVERAGE ABUNDANCE\n')
fprintf('sec\tNadjSPC\tCOV\n')
for j=1:4
    ind = find(strcmp(ptd,'yes') & NadjSPC(:,j)>=cutoff);
    fprintf('%s\t%f\t%f\n',sec{j},mean(NadjSPC(ind,j)),mean(COV(ind,j)));
end
% -- non-plastid genes --
fprintf('\n-- Non-Plastid proteins --\n');
fprintf('\nCORRELATION USING LOG10\n')
fprintf('type\t#prot\tRp\tp.val\tRs\tp.val\n')
ind_np = find(~strcmp(ptd,'yes'));
[Rp,Pp] = corr(log10(mean(NadjSPC(ind_np,:),2)), log10(mean(COV(ind_np,:),2)), 'type', 'Pearson');
[Rs,Ps] = corr(log10(mean(NadjSPC(ind_np,:),2)), log10(mean(COV(ind_np,:),2)), 'type', 'Spearman');
fprintf('average\t%d\t%f\t%f\t%f\t%f\n',length(ind_np),Rp,Pp,Rs,Ps);
for j=1:4
    ind = find(~strcmp(ptd,'yes') & NadjSPC(:,j)>=cutoff);
    [rp,pp] = corr(log10(NadjSPC(ind,j)),log10(COV(ind,j)),'type','Pearson');
    [rs,ps] = corr(log10(NadjSPC(ind,j)),log10(COV(ind,j)),'type','Spearman');
    fprintf('%s\t%d\t%f\t%f\t%f\t%f\n',sec{j},length(ind),rp,pp,rs,ps);
end
fprintf('\nAVERAGE ABUNDANCE\n')
fprintf('sec\tNadjSPC\tCOV\n')
for j=1:4
    ind = find(~strcmp(ptd,'yes') & NadjSPC(:,j)>=cutoff);
    fprintf('%s\t%f\t%f\n',sec{j},mean(NadjSPC(ind,j)),mean(COV(ind,j)));
end


fprintf('\n NSAF & RPKM \n')
% -- plastid genes --
fprintf('\n-- Plastid proteins --\n');
fprintf('\nCORRELATION USING LOG10\n')
fprintf('type\t#prot\tRp\tp.val\tRs\tp.val\n')
ind_p = find(strcmp(ptd,'yes'));
[Rp,Pp] = corr(log10(mean(NSAF(ind_p,:),2)), log10(mean(RPKM(ind_p,:),2)), 'type', 'Pearson');
[Rs,Ps] = corr(log10(mean(NSAF(ind_p,:),2)), log10(mean(RPKM(ind_p,:),2)), 'type', 'Spearman');
fprintf('average\t%d\t%f\t%f\t%f\t%f\n',length(ind_p),Rp,Pp,Rs,Ps);
for j=1:4
    ind = find(strcmp(ptd,'yes') & NadjSPC(:,j)>=cutoff);
    [rp,pp] = corr(log10(NSAF(ind,j)),log10(RPKM(ind,j)),'type','Pearson');
    [rs,ps] = corr(log10(NSAF(ind,j)),log10(RPKM(ind,j)),'type','Spearman');
    fprintf('%s\t%d\t%f\t%f\t%f\t%f\n',sec{j},length(ind),rp,pp,rs,ps);
end
fprintf('\nAVERAGE ABUNDANCE\n')
fprintf('sec\tNSAF\tRPKM\n')
for j=1:4
    ind = find(strcmp(ptd,'yes') & NadjSPC(:,j)>=cutoff);
    fprintf('%s\t%f\t%f\n',sec{j},mean(NSAF(ind,j)),mean(RPKM(ind,j)));
end
% -- non-plastid genes --
fprintf('\n-- Non-Plastid proteins --\n');
fprintf('\nCORRELATION USING LOG10\n')
fprintf('type\t#prot\tRp\tp.val\tRs\tp.val\n')
ind_np = find(~strcmp(ptd,'yes'));
[Rp,Pp] = corr(log10(mean(NSAF(ind_np,:),2)), log10(mean(RPKM(ind_np,:),2)), 'type', 'Pearson');
[Rs,Ps] = corr(log10(mean(NSAF(ind_np,:),2)), log10(mean(RPKM(ind_np,:),2)), 'type', 'Spearman');
fprintf('average\t%d\t%f\t%f\t%f\t%f\n',length(ind_np),Rp,Pp,Rs,Ps);
for j=1:4
    ind = find(~strcmp(ptd,'yes') & NadjSPC(:,j)>=cutoff);
    [rp,pp] = corr(log10(NSAF(ind,j)),log10(RPKM(ind,j)),'type','Pearson');
    [rs,ps] = corr(log10(NSAF(ind,j)),log10(RPKM(ind,j)),'type','Spearman');
    fprintf('%s\t%d\t%f\t%f\t%f\t%f\n',sec{j},length(ind),rp,pp,rs,ps);
end
fprintf('\nAVERAGE ABUNDANCE\n')
fprintf('sec\tNSAF\tRPKM\n')
for j=1:4
    ind = find(~strcmp(ptd,'yes') & NadjSPC(:,j)>=cutoff);
    fprintf('%s\t%f\t%f\n',sec{j},mean(NSAF(ind,j)),mean(RPKM(ind,j)));
end

% % -- generate plots --
% % save('c1p.mat','c1p'); save('c1np.mat','c1np'); save('c2p.mat','c2p'); save('c2np.mat','c2np')
% load 'c1p.mat'; load 'c1np.mat'; load 'c2p.mat'; load 'c2np.mat'
% figure, plot(c1p(:,2),c1p(:,3),'ob','MarkerFaceColor','b'), hold on, plot(c1np(:,2),c1np(:,3),'or','MarkerFaceColor','r'), hold off,  grid, legend('plastid','non-plastid','Location','NorthWest'), xlabel('NadjSPC'), ylabel('COV'), title('Average abundance'), text(c1p(:,2),c1p(:,3),sec), text(c1np(:,2),c1np(:,3),sec), saveas(gcf,'cut_NadjSPC-COV','jpg')
% figure, plot(c2p(:,2),c2p(:,3),'ob','MarkerFaceColor','b'), hold on, plot(c2np(:,2),c2np(:,3),'or','MarkerFaceColor','r'), hold off, grid, legend('plastid','non-plastid','Location','NorthWest'), xlabel('NSAF'), ylabel('RPKM'), title('Average abundance'), text(c2p(:,2),c2p(:,3),sec), text(c2np(:,2),c2np(:,3),sec), saveas(gcf,'cut_NSAF-RPKM','jpg')
% 
% % save('d1p.mat','d1p'); save('d1np.mat','d1np'); save('d2p.mat','d2p'); save('d2np.mat','d2np')
% load 'd1p.mat'; load 'd1np.mat'; load 'd2p.mat'; load 'd2np.mat'
% figure, plot(d1p(:,2),d1p(:,3),'ob','MarkerFaceColor','b'), hold on, plot(d1np(:,2),d1np(:,3),'or','MarkerFaceColor','r'), hold off, grid, legend('plastid','non-plastid','Location','NorthWest'), xlabel('NadjSPC'), ylabel('COV'), title('Average abundance'), text(d1p(:,2),d1p(:,3),sec), text(d1np(:,2),d1np(:,3),sec), saveas(gcf,'all_NadjSPC-COV','jpg')
% figure, plot(d2p(:,2),d2p(:,3),'ob','MarkerFaceColor','b'), hold on, plot(d2np(:,2),d2np(:,3),'or','MarkerFaceColor','r'), hold off, grid, legend('plastid','non-plastid','Location','NorthWest'), xlabel('NSAF'), ylabel('RPKM'), title('Average abundance'), text(d2p(:,2),d2p(:,3),sec), text(d2np(:,2),d2np(:,3),sec), saveas(gcf,'all_NSAF-RPKM','jpg')
