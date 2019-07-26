clear;

% -----------------------------------------
% Fig2 (revised version)
% extract outliers as marked
% -----------------------------------------

[N,T] = xlsread('Supplemental Table 1A-B-C-D-E-aug2013-v1.xls',5);
NSAF = N(:,4:7);
RPKM = N(:,8:11);
NadjSPC = N(:,13:16);
COV = N(:,17:20);
gid = cellstr(num2str(N(:,1))); % T(3:end,2);

cutoff = 0.00005;
colors = ['r','g','b','k'];

tklabSize=22; tklabWeight='normal'; tklabFont='Arial';
axlabSize=22; axlabWeight='normal'; axlabFont='Arial';


% fprintf('\n NadjSPC & COV \n')

% average
NadjSPC_avg = mean(NadjSPC,2); COV_avg = mean(COV,2);
ind = find(NadjSPC_avg>=cutoff);
figure, plot(log10(NadjSPC_avg(ind)),log10(COV_avg(ind)),'b*'),
text(log10(NadjSPC_avg(ind)),log10(COV_avg(ind)),gid(ind),'horizontal','left', 'vertical','bottom')
fprintf('Number of proteins = %d\n',length(ind));
axis([-4.5,-1.5,0,6])
xlabel('log10(NadjSPC)'); ylabel('log10(COV)');
% set(gca,'FontSize',tklabSize,'FontWeight',tklabWeight,'FontName',tklabFont)
% set(findall(gcf,'type','text'),'FontSize',axlabSize,'FontWeight',axlabWeight,'FontName',axlabFont)
saveas(gcf, ['Fig2A.avg_NadjSPC_',num2str(cutoff),'.jpg'], 'jpg');

% section-wise
figure
for j=1:4
    ind = find(NadjSPC(:,j)>=cutoff);
    fprintf('sec %d : %d proteins\n',j,length(ind));
    plot(log10(NadjSPC(ind,j)),log10(COV(ind,j)),[colors(j),'.'])
    text(log10(NadjSPC(ind,j)),log10(COV(ind,j)),gid(ind),'horizontal','left', 'vertical','bottom')
    hold on
end
hold off
axis([-4.5,-1.5,0,6])
xlabel('log10(NadjSPC)'); ylabel('log10(COV)');
legend('0-1cm','3-4cm','8-9cm','13-14cm','Location','NorthWest')
% set(gca,'FontSize',tklabSize,'FontWeight',tklabWeight,'FontName',tklabFont)
% set(findall(gcf,'type','text'),'FontSize',axlabSize,'FontWeight',axlabWeight,'FontName',axlabFont)
saveas(gcf,['Fig2C.sec_NadjSPC_',num2str(cutoff),'.jpg'],'jpg')


% fprintf('\n NSAF & RPKM \n')

% average
NSAF_avg = mean(NSAF,2); RPKM_avg = mean(RPKM,2);
ind = find(NadjSPC_avg>=cutoff);
figure, plot(log10(NSAF_avg(ind)),log10(RPKM_avg(ind)),'b*'),
text(log10(NSAF_avg(ind)),log10(RPKM_avg(ind)),gid(ind),'horizontal','left', 'vertical','bottom')
fprintf('Number of proteins = %d\n',length(ind));
axis([-5.5,-1.5,-1,5])
xlabel('log10(NSAF)'); ylabel('log10(RPKM)');
% set(gca,'FontSize',tklabSize,'FontWeight',tklabWeight,'FontName',tklabFont)
% set(findall(gcf,'type','text'),'FontSize',axlabSize,'FontWeight',axlabWeight,'FontName',axlabFont)
saveas(gcf, ['Fig2B.avg_NSAF_',num2str(cutoff),'.jpg'], 'jpg');

% section-wise
figure
for j=1:4
    ind = find(NadjSPC(:,j)>=cutoff);
    fprintf('sec %d : %d proteins\n',j,length(ind));
    plot(log10(NSAF(ind,j)),log10(RPKM(ind,j)),[colors(j),'.'])
    text(log10(NSAF(ind,j)),log10(RPKM(ind,j)),gid(ind),'horizontal','left', 'vertical','bottom')
    hold on
end
hold off
axis([-5.5,-1.5,-1,5])
xlabel('log10(NSAF)'); ylabel('log10(RPKM)');
legend('0-1cm','3-4cm','8-9cm','13-14cm','Location','NorthWest')
% set(gca,'FontSize',tklabSize,'FontWeight',tklabWeight,'FontName',tklabFont)
% set(findall(gcf,'type','text'),'FontSize',axlabSize,'FontWeight',axlabWeight,'FontName',axlabFont)
saveas(gcf,['Fig2D.sec_NSAF_',num2str(cutoff),'.jpg'],'jpg')

% -----------------------------------------
