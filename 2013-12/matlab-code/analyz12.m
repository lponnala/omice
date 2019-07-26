clear;

[N,T] = xlsread('masterclusterv12-80Slalit.xlsx');
NSAF = N(:,19:22);
RPKM = N(:,23:26);

sec = {'sec1','sec4','sec9','sec14'};
symb = {'bd','rs','g^','ko'};

nsaf_all = []; rpkm_all = [];
for j=1:4
    ind = find(NSAF(:,j)>0);
    nsaf_all = [nsaf_all; NSAF(ind,j)]; rpkm_all = [rpkm_all; RPKM(ind,j)];
    plot(log10(NSAF(ind,j)),log10(RPKM(ind,j)),symb{j});
    hold on
end
hold off
legend('0-1cm','3-4cm','9-10cm','13-14cm');

[Rp,Pp] = corr(log10(nsaf_all),log10(rpkm_all),'type','Pearson');
[Rs,Ps] = corr(log10(nsaf_all),log10(rpkm_all),'type','Spearman');
fprintf('number of proteins = %d\n',length(nsaf_all));
fprintf('Pearson:\ncorr = %f, pvalue = %f\n',Rp,Pp);
fprintf('Spearman:\ncorr = %f, pvalue = %f\n',Rs,Ps);
