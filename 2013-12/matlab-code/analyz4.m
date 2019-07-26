% // Investigate non-linearity //
clear;

[N,T] = xlsread('Supplemental Table 1A-B-C-D-E-aug2013-v1.xls',5);
NSAF = N(:,4:7);
RPKM = N(:,8:11);
NadjSPC = N(:,13:16);
COV = N(:,17:20);

j=4; % section (1,2,3,4)
nadjspc = NadjSPC(:,j);
cov = COV(:,j);

figure
subplot(221), plot(log10(nadjspc),log10(cov),'b.')
title('All proteins'), xlabel('log10(NadjSPC)'), ylabel('log10(COV)')
ind = find(nadjspc>=quantile(nadjspc,0.9));
subplot(222), plot(log10(nadjspc(ind)),log10(cov(ind)),'b.')
title(['Top: ',num2str(length(ind)),' proteins']), xlabel('log10(NadjSPC)'), ylabel('log10(COV)')
ind = find(nadjspc>=quantile(nadjspc,0.7) & nadjspc<quantile(nadjspc,0.8));
subplot(223), plot(log10(nadjspc(ind)),log10(cov(ind)),'b.')
title(['Middle: ',num2str(length(ind)),' proteins']), xlabel('log10(NadjSPC)'), ylabel('log10(COV)')
ind = find(nadjspc>=quantile(nadjspc,0.4) & nadjspc<quantile(nadjspc,0.5));
subplot(224), plot(log10(nadjspc(ind)),log10(cov(ind)),'b.')
title(['Bottom: ',num2str(length(ind)),' proteins']), xlabel('log10(NadjSPC)'), ylabel('log10(COV)')
saveas(gcf,['sec',num2str(j)],'jpg')

