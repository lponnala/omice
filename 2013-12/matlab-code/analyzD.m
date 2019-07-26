clear;

% ----------------------------------------------------------------
% Using all GRMZM genes
% ----------------------------------------------------------------
[N,T] = xlsread('mRNAabundance-freqObservation-lalit.xlsx');

COV=N(:,1); RPKM=N(:,2);
NSAF=N(:,3); NadjSPC=N(:,4);
gid = T(2:end,2);

detect_thresh = 0;
num_bins = 100;

x = RPKM; y = NSAF; x_label='RPKM'; y_label='fraction of detected proteins (NSAF>0)';

qt = quantile(x,num_bins);
det_ratio = nan(1,num_bins);
ind = find(x<=qt(1));
det_ratio(1) = length(find(y(ind)>detect_thresh))/length(ind); % ignores NaN values
fprintf('1\t%.4f\t%d\t%d\t%.4f\n',qt(1),length(find(y(ind)>detect_thresh)),length(ind),det_ratio(1));
for k=2:num_bins
    ind = find(x>qt(k-1) & x<=qt(k));
    det_ratio(k) = length(find(y(ind)>detect_thresh))/length(ind);
    fprintf('%d\t%.2f\t%d\t%d\t%.4f\n',k,qt(k),length(find(y(ind)>detect_thresh)),length(ind),det_ratio(k));
end
plot(qt,det_ratio), hold on
hold off
xlabel(x_label); ylabel(y_label);
% ----------------------------------------------------------------


% % ----------------------------------------------------------------
% % Using matched proteins
% % ----------------------------------------------------------------
% 
% [N,T] = xlsread('forLalit-SupplTable1D-correlations.xlsx');
% NSAF = N(:,5:8);
% NadjSPC = N(:,9:12);
% COV = N(:,14:17);
% RPKM = N(:,18:21);
% 
% % [N,T] = xlsread('Supplemental Table 1A-B-C-D-E-aug2013-v1.xls',5);
% % NSAF = N(:,4:7);
% % RPKM = N(:,8:11);
% % NadjSPC = N(:,13:16);
% % COV = N(:,17:20);
% 
% sec = {'sec1','sec4','sec9','sec14'};
% detect_thresh = 0;
% num_bins = 10;
% 
% % X = COV; Y = NadjSPC; X_label='COV'; Y_label='fraction of detected proteins (NadjSPC>0)';
% X = RPKM; Y = NSAF; X_label='RPKM'; Y_label='fraction of detected proteins (NSAF>0)';
% 
% color = {'b','r','g','k'};
% figure
% for j=1:4
%     fprintf('\n-- %s --\n',sec{j});
%     ind = find(X(:,j)>0);
%     x = X(ind,j); y = Y(ind,j);
%     qt = quantile(x,num_bins);
%     det_ratio = nan(1,num_bins);
%     ind = find(x<=qt(1));
%     det_ratio(1) = length(find(y(ind)>detect_thresh))/length(ind);
%     fprintf('1\t%.2f\t%d\t%d\t%.4f\n',qt(1),length(find(y(ind)>detect_thresh)),length(ind),det_ratio(1));
%     for k=2:num_bins
%         ind = find(x>qt(k-1) & x<=qt(k));
%         det_ratio(k) = length(find(y(ind)>detect_thresh))/length(ind);
%         fprintf('%d\t%.2f\t%d\t%d\t%.4f\n',k,qt(k),length(find(y(ind)>detect_thresh)),length(ind),det_ratio(k));
%     end
%     plot(qt,det_ratio,color{j}), hold on
% end
% hold off
% legend(sec,'Location','Best')
% xlabel(X_label); ylabel(Y_label);
% 
% 
% % for j=1:4
% %     ind = find(X(:,j)>0);
% %     x = X(ind,j); y = Y(ind,j);
% %     qt = quantile(x,num_bins);
% %     det_ratio = nan(1,num_bins);
% %     ind = find(x<=qt(1));
% %     det_ratio(1) = length(find(y(ind)>detect_thresh))/length(ind);
% %     for k=2:num_bins
% %         ind = find(x>qt(k-1) & x<=qt(k));
% %         det_ratio(k) = length(find(y(ind)>detect_thresh))/length(ind);
% %     end
% %     figure, plot(qt,det_ratio)
% %     title(sec{j})
% % end
% 
% % ----------------------------------------------------------------

