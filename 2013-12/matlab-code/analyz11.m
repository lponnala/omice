clear;

% % -----------------------------------------
% % Fig1 (revised version)
% % -----------------------------------------
% 
% [N,T]=xlsread('newsuppltable1A.xlsx');
% cov=N(:,5:8); rpkm=N(:,9:12);
% type=T(2:end,4); match=T(2:end,1);
% geneid = T(2:end,6); geneidT = T(2:end,7); geneidP = T(2:end,8);
% 
% colors = ['r','g','b','k'];
% 
% tklabSize=22; tklabWeight='normal'; tklabFont='Arial';
% axlabSize=22; axlabWeight='normal'; axlabFont='Arial';
% 
% ind = find(strncmp('GRMZM',geneid,5));
% fprintf('A,B : number of genes = %d\n',length(ind));
% COV=cov(ind,:); RPKM=rpkm(ind,:);
% figure;
% for j=1:4
%     ind=find(RPKM(:,j)>0);
%     data=log10(RPKM(ind,j));
%     [nelements,centers] = hist(data);
%     plot(centers,nelements/sum(nelements),['-o',colors(j)],'MarkerFaceColor',colors(j));
%     hold on
% end
% xlabel('log10(RPKM)'), ylabel('Frequency')
% axis([-3,4,0,0.4])
% legend('0-1cm','3-4cm','8-9cm','13-14cm','Location','NorthEast')
% hold off
% set(gca,'FontSize',tklabSize,'FontWeight',tklabWeight,'FontName',tklabFont)
% set(findall(gcf,'type','text'),'FontSize',axlabSize,'FontWeight',axlabWeight,'FontName',axlabFont)
% saveas(gcf,'Fig1A','jpg')
% figure;
% for j=1:4
%     ind=find(COV(:,j)>0);
%     data=log10(COV(ind,j));
%     [nelements,centers] = hist(data);
%     plot(centers,nelements/sum(nelements),['-o',colors(j)],'MarkerFaceColor',colors(j));
%     hold on
% end
% xlabel('log10(COV)'), ylabel('Frequency')
% axis([-1,6,0,0.4])
% legend('0-1cm','3-4cm','8-9cm','13-14cm','Location','NorthEast')
% hold off
% set(gca,'FontSize',tklabSize,'FontWeight',tklabWeight,'FontName',tklabFont)
% set(findall(gcf,'type','text'),'FontSize',axlabSize,'FontWeight',axlabWeight,'FontName',axlabFont)
% saveas(gcf,'Fig1B','jpg')
% 
% ind = find(strncmp('GRMZM',geneid,5) & strncmp('f',type,1));
% fprintf('C,D : number of genes = %d\n',length(ind));
% COV=cov(ind,:); RPKM=rpkm(ind,:);
% figure;
% for j=1:4
%     ind=find(RPKM(:,j)>0);
%     data=log10(RPKM(ind,j));
%     [nelements,centers] = hist(data);
%     plot(centers,nelements/sum(nelements),['-o',colors(j)],'MarkerFaceColor',colors(j));
%     hold on
% end
% xlabel('log10(RPKM)'), ylabel('Frequency')
% axis([-3,4,0,0.4])
% legend('0-1cm','3-4cm','8-9cm','13-14cm','Location','NorthEast')
% hold off
% set(gca,'FontSize',tklabSize,'FontWeight',tklabWeight,'FontName',tklabFont)
% set(findall(gcf,'type','text'),'FontSize',axlabSize,'FontWeight',axlabWeight,'FontName',axlabFont)
% saveas(gcf,'Fig1C','jpg')
% figure;
% for j=1:4
%     ind=find(COV(:,j)>0);
%     data=log10(COV(ind,j));
%     [nelements,centers] = hist(data);
%     plot(centers,nelements/sum(nelements),['-o',colors(j)],'MarkerFaceColor',colors(j));
%     hold on
% end
% xlabel('log10(COV)'), ylabel('Frequency')
% axis([-1,6,0,0.4])
% legend('0-1cm','3-4cm','8-9cm','13-14cm','Location','NorthEast')
% hold off
% set(gca,'FontSize',tklabSize,'FontWeight',tklabWeight,'FontName',tklabFont)
% set(findall(gcf,'type','text'),'FontSize',axlabSize,'FontWeight',axlabWeight,'FontName',axlabFont)
% saveas(gcf,'Fig1D','jpg')
% 
% ind = find(strncmp('GRMZM',geneid,5) & strncmp('x',match,1));
% fprintf('E,F : number of genes = %d\n',length(ind));
% COV=cov(ind,:); RPKM=rpkm(ind,:);
% figure;
% for j=1:4
%     ind=find(RPKM(:,j)>0);
%     data=log10(RPKM(ind,j));
%     [nelements,centers] = hist(data);
%     plot(centers,nelements/sum(nelements),['-o',colors(j)],'MarkerFaceColor',colors(j));
%     hold on
% end
% xlabel('log10(RPKM)'), ylabel('Frequency')
% axis([-3,4,0,0.4])
% legend('0-1cm','3-4cm','8-9cm','13-14cm','Location','NorthEast')
% hold off
% set(gca,'FontSize',tklabSize,'FontWeight',tklabWeight,'FontName',tklabFont)
% set(findall(gcf,'type','text'),'FontSize',axlabSize,'FontWeight',axlabWeight,'FontName',axlabFont)
% saveas(gcf,'Fig1E','jpg')
% figure;
% for j=1:4
%     ind=find(COV(:,j)>0);
%     data=log10(COV(ind,j));
%     [nelements,centers] = hist(data);
%     plot(centers,nelements/sum(nelements),['-o',colors(j)],'MarkerFaceColor',colors(j));
%     hold on
% end
% xlabel('log10(COV)'), ylabel('Frequency')
% axis([-1,6,0,0.4])
% legend('0-1cm','3-4cm','8-9cm','13-14cm','Location','NorthEast')
% hold off
% set(gca,'FontSize',tklabSize,'FontWeight',tklabWeight,'FontName',tklabFont)
% set(findall(gcf,'type','text'),'FontSize',axlabSize,'FontWeight',axlabWeight,'FontName',axlabFont)
% saveas(gcf,'Fig1F','jpg')
% 
% [N,T] = xlsread('Supplemental Table 1A-B-C-D-E-aug2013-v1.xls',5);
% RPKM = N(:,8:11); COV = N(:,17:20);
% fprintf('G,H : Number of genes = %d\n',size(RPKM,1));
% figure;
% for j=1:4
%     ind=find(RPKM(:,j)>0);
%     data=log10(RPKM(ind,j));
%     [nelements,centers] = hist(data);
%     plot(centers,nelements/sum(nelements),['-o',colors(j)],'MarkerFaceColor',colors(j));
%     hold on
% end
% xlabel('log10(RPKM)'), ylabel('Frequency')
% axis([-3,4,0,0.4])
% legend('0-1cm','3-4cm','8-9cm','13-14cm','Location','NorthWest')
% hold off
% set(gca,'FontSize',tklabSize,'FontWeight',tklabWeight,'FontName',tklabFont)
% set(findall(gcf,'type','text'),'FontSize',axlabSize,'FontWeight',axlabWeight,'FontName',axlabFont)
% saveas(gcf,'Fig1G','jpg')
% figure;
% for j=1:4
%     ind=find(COV(:,j)>0);
%     data=log10(COV(ind,j));
%     [nelements,centers] = hist(data);
%     plot(centers,nelements/sum(nelements),['-o',colors(j)],'MarkerFaceColor',colors(j));
%     hold on
% end
% xlabel('log10(COV)'), ylabel('Frequency')
% axis([-1,6,0,0.4])
% legend('0-1cm','3-4cm','8-9cm','13-14cm','Location','NorthEast')
% hold off
% set(gca,'FontSize',tklabSize,'FontWeight',tklabWeight,'FontName',tklabFont)
% set(findall(gcf,'type','text'),'FontSize',axlabSize,'FontWeight',axlabWeight,'FontName',axlabFont)
% saveas(gcf,'Fig1H','jpg')
% 
% [N,T]=xlsread('DataForFig1GH.xlsx');
% NSAF=N(:,7:10); NadjSPC=N(:,15:18);
% fprintf('I,J : Number of genes = %d\n',size(NSAF,1));
% figure;
% for j=1:4
%     ind=find(NSAF(:,j)>0);
%     data=log10(NSAF(ind,j));
%     [nelements,centers] = hist(data);
%     plot(centers,nelements/sum(nelements),['-o',colors(j)],'MarkerFaceColor',colors(j));
%     hold on
% end
% xlabel('log10(NSAF)'), ylabel('Frequency')
% axis([-6,-1,0,0.25])
% legend('0-1cm','3-4cm','8-9cm','13-14cm','Location','NorthWest')
% hold off
% set(gca,'FontSize',tklabSize,'FontWeight',tklabWeight,'FontName',tklabFont)
% set(findall(gcf,'type','text'),'FontSize',axlabSize,'FontWeight',axlabWeight,'FontName',axlabFont)
% saveas(gcf,'Fig1I','jpg')
% figure;
% for j=1:4
%     ind=find(NadjSPC(:,j)>0);
%     data=log10(NadjSPC(ind,j));
%     [nelements,centers] = hist(data);
%     plot(centers,nelements/sum(nelements),['-o',colors(j)],'MarkerFaceColor',colors(j));
%     hold on
% end
% xlabel('log10(NadjSPC)'), ylabel('Frequency')
% axis([-6,-1,0,0.35])
% legend('0-1cm','3-4cm','8-9cm','13-14cm','Location','NorthEast')
% hold off
% set(gca,'FontSize',tklabSize,'FontWeight',tklabWeight,'FontName',tklabFont)
% set(findall(gcf,'type','text'),'FontSize',axlabSize,'FontWeight',axlabWeight,'FontName',axlabFont)
% saveas(gcf,'Fig1J','jpg')
% 
% [N,T] = xlsread('Supplemental Table 1A-B-C-D-E-aug2013-v1.xls',5);
% NSAF = N(:,4:7); NadjSPC = N(:,13:16);
% fprintf('K,L : Number of genes = %d\n',size(NSAF,1));
% figure;
% for j=1:4
%     ind=find(NSAF(:,j)>0);
%     data=log10(NSAF(ind,j));
%     [nelements,centers] = hist(data);
%     plot(centers,nelements/sum(nelements),['-o',colors(j)],'MarkerFaceColor',colors(j));
%     hold on
% end
% xlabel('log10(NSAF)'), ylabel('Frequency')
% axis([-6,-1,0,0.25])
% legend('0-1cm','3-4cm','8-9cm','13-14cm','Location','NorthWest')
% hold off
% set(gca,'FontSize',tklabSize,'FontWeight',tklabWeight,'FontName',tklabFont)
% set(findall(gcf,'type','text'),'FontSize',axlabSize,'FontWeight',axlabWeight,'FontName',axlabFont)
% saveas(gcf,'Fig1K','jpg')
% figure;
% for j=1:4
%     ind=find(NadjSPC(:,j)>0);
%     data=log10(NadjSPC(ind,j));
%     [nelements,centers] = hist(data);
%     plot(centers,nelements/sum(nelements),['-o',colors(j)],'MarkerFaceColor',colors(j));
%     hold on
% end
% xlabel('log10(NadjSPC)'), ylabel('Frequency')
% axis([-6,-1,0,0.35])
% legend('0-1cm','3-4cm','8-9cm','13-14cm','Location','NorthEast')
% hold off
% set(gca,'FontSize',tklabSize,'FontWeight',tklabWeight,'FontName',tklabFont)
% set(findall(gcf,'type','text'),'FontSize',axlabSize,'FontWeight',axlabWeight,'FontName',axlabFont)
% saveas(gcf,'Fig1L','jpg')
% 
% close all
% 
% % [N,T]=xlsread('suppltable1Sept10.xlsx');
% % COV=N(:,2:5); RPKM=N(:,6:9);
% % type=T(2:end,4); match=T(2:end,1);
% % geneid = T(2:end,6); geneidT = T(2:end,7); geneidP = T(2:end,8);
% 
% % % // verifying numbers in the current Fig.1 //
% % % A,B: sum(strncmp('GRMZM',geneid,5)) = 33498
% % % C,D: sum(strncmp('GRMZM',geneid,5) & strncmp('x',match,1)) = 2813
% % % E,F: [N,T] = xlsread('Supplemental Table 1A-B-C-D-E-aug2013-v1.xls',5); RPKM = N(:,8:11); COV = N(:,17:20);
% % % G,H: ??
% % % I,J: [N,T] = xlsread('Supplemental Table 1A-B-C-D-E-aug2013-v1.xls',5); NSAF = N(:,4:7); NadjSPC = N(:,13:16);
% % 
% % % // new Fig.1 //
% % % A,B: sum(strncmp('GRMZM',geneid,5) & strncmp('f',type,1)) = 22392
% 
% % -----------------------------------------


% -----------------------------------------
% Fig2 (revised version)
% -----------------------------------------

[N,T] = xlsread('Supplemental Table 1A-B-C-D-E-aug2013-v1.xls',5);
NSAF = N(:,4:7);
RPKM = N(:,8:11);
NadjSPC = N(:,13:16);
COV = N(:,17:20);

cutoff = 0.00005;
colors = ['r','g','b','k'];

tklabSize=22; tklabWeight='normal'; tklabFont='Arial';
axlabSize=22; axlabWeight='normal'; axlabFont='Arial';

fprintf('\n NadjSPC & COV \n')
% average
NadjSPC_avg = mean(NadjSPC,2); COV_avg = mean(COV,2);
ind = find(NadjSPC_avg>=cutoff);
figure, plot(log10(NadjSPC_avg(ind)),log10(COV_avg(ind)),'b*');
fprintf('Number of proteins = %d\n',length(ind));
axis([-4.5,-1.5,0,6])
xlabel('log10(NadjSPC)'); ylabel('log10(COV)');
set(gca,'FontSize',tklabSize,'FontWeight',tklabWeight,'FontName',tklabFont)
set(findall(gcf,'type','text'),'FontSize',axlabSize,'FontWeight',axlabWeight,'FontName',axlabFont)
saveas(gcf, ['Fig2A.avg_NadjSPC_',num2str(cutoff),'.jpg'], 'jpg');
% section-wise
figure
for j=1:4
    ind = find(NadjSPC(:,j)>=cutoff);
    fprintf('sec %d : %d proteins\n',j,length(ind));
    plot(log10(NadjSPC(ind,j)),log10(COV(ind,j)),[colors(j),'.'])
    hold on
end
hold off
axis([-4.5,-1.5,0,6])
xlabel('log10(NadjSPC)'); ylabel('log10(COV)');
legend('0-1cm','3-4cm','8-9cm','13-14cm','Location','NorthWest')
set(gca,'FontSize',tklabSize,'FontWeight',tklabWeight,'FontName',tklabFont)
set(findall(gcf,'type','text'),'FontSize',axlabSize,'FontWeight',axlabWeight,'FontName',axlabFont)
saveas(gcf,['Fig2C.sec_NadjSPC_',num2str(cutoff),'.jpg'],'jpg')

fprintf('\n NSAF & RPKM \n')
% average
NSAF_avg = mean(NSAF,2); RPKM_avg = mean(RPKM,2);
ind = find(NadjSPC_avg>=cutoff);
figure, plot(log10(NSAF_avg(ind)),log10(RPKM_avg(ind)),'b*');
fprintf('Number of proteins = %d\n',length(ind));
axis([-5.5,-1.5,-1,5])
xlabel('log10(NSAF)'); ylabel('log10(RPKM)');
set(gca,'FontSize',tklabSize,'FontWeight',tklabWeight,'FontName',tklabFont)
set(findall(gcf,'type','text'),'FontSize',axlabSize,'FontWeight',axlabWeight,'FontName',axlabFont)
saveas(gcf, ['Fig2B.avg_NSAF_',num2str(cutoff),'.jpg'], 'jpg');
% section-wise
figure
for j=1:4
    ind = find(NadjSPC(:,j)>=cutoff);
    fprintf('sec %d : %d proteins\n',j,length(ind));
    plot(log10(NSAF(ind,j)),log10(RPKM(ind,j)),[colors(j),'.'])
    hold on
end
hold off
axis([-5.5,-1.5,-1,5])
xlabel('log10(NSAF)'); ylabel('log10(RPKM)');
legend('0-1cm','3-4cm','8-9cm','13-14cm','Location','NorthWest')
set(gca,'FontSize',tklabSize,'FontWeight',tklabWeight,'FontName',tklabFont)
set(findall(gcf,'type','text'),'FontSize',axlabSize,'FontWeight',axlabWeight,'FontName',axlabFont)
saveas(gcf,['Fig2D.sec_NSAF_',num2str(cutoff),'.jpg'],'jpg')

% -----------------------------------------


% % -----------------------------------------
% % Fig 2U (new)
% % -----------------------------------------
% 
% [N,T]=xlsread('newsuppltable1A.xlsx');
% cov=N(:,5:8); rpkm=N(:,9:12);
% type=T(2:end,4); match=T(2:end,1);
% geneid = T(2:end,6); geneidT = T(2:end,7); geneidP = T(2:end,8);
% 
% colors = ['r','g','b','k'];
% ind = find(strncmp('GRMZM',geneid,5) & strncmp('x',match,1));
% fprintf('C,D : number of genes = %d\n',length(ind));
% COV=cov(ind,:); RPKM=rpkm(ind,:);
% figure;
% for j=1:4
%     ind=find(RPKM(:,j)>0);
%     data=log10(RPKM(ind,j));
%     [nelements,centers] = hist(data);
%     plot(centers,nelements/sum(nelements),['-o',colors(j)],'MarkerFaceColor',colors(j));
%     hold on
% end
% xlabel('log10(RPKM)'), ylabel('Frequency')
% legend('0-1cm','3-4cm','8-9cm','13-14cm','Location','NorthEast')
% hold off
% saveas(gcf,'Fig1C','jpg')
% figure;
% for j=1:4
%     ind=find(COV(:,j)>0);
%     data=log10(COV(ind,j));
%     [nelements,centers] = hist(data);
%     plot(centers,nelements/sum(nelements),['-o',colors(j)],'MarkerFaceColor',colors(j));
%     hold on
% end
% xlabel('log10(COV)'), ylabel('Frequency')
% legend('0-1cm','3-4cm','8-9cm','13-14cm','Location','NorthEast')
% hold off
% saveas(gcf,'Fig1D','jpg')
% 
% % -----------------------------------------

