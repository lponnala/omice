% // Clustering protein/mRNA ratio at the gene level //
clear;

% ----------------------------------------------------
% Check for unique assignment of genes to clusters
% ----------------------------------------------------
file = 'NadjSPC_COV_repsmall.xlsx'; % 'NSAF_RPKM.xlsx', 'NSAF_RPKM_repsmall.xlsx', 'NadjSPC_COV.xlsx', 'NadjSPC_COV_repsmall.xlsx'
[N,T]=xlsread(file);
ids = N(2:end,:);
for j=1:(size(ids,2)-1)
    for k=(j+1):size(ids,2)
        fprintf('-- comparing clusters %d, %d --\n',j,k);
        % check if any of the ids in column j are present in column k
        indj = find(isfinite(ids(:,j))); idj = ids(indj,j);
        indk = find(isfinite(ids(:,k))); idk = ids(indk,k);
        fprintf('length(idj) = %d, length(idk) = %d\n',length(idj),length(idk));
        for i=1:length(idj)
            if ~isempty(find(idk==idj(i)))
                error('matching id %d found in columns %d, %d\n',idj(i),j,k); break
            end
        end
    end
end
% ----------------------------------------------------


% % ----------------------------------------------------
% % Create the dendrogram
% % ----------------------------------------------------
% [N,T] = xlsread('Supplemental Table 1A-B-C-D-E-aug2013-v1.xls',5);
% NSAF = N(:,4:7);
% RPKM = N(:,8:11);
% NadjSPC = N(:,13:16);
% COV = N(:,17:20);
% groupID = N(:,1);
% 
% num = NSAF; % NadjSPC;
% den = RPKM; % COV;
% 
% % % replace 0 values of num with 1e-6
% % num(find(num==0))=1e-6;
% 
% % filter out rows containing NaN or Inf
% ratio = num./den;
% ind = find(isfinite(sum(ratio,2)));
% ratio = ratio(ind,:);
% groupID = groupID(ind,:);
% cols = {'sec1','sec4','sec9','sec14'};
% 
% outname='nsaf_ratio';
% CGobj=clustergram(ratio,'RowLabels',groupID,'ColumnLabels',cols,'RowPdist','correlation','Linkage','average','Cluster',1);
% plot(CGobj); title([outname,' clusters']); saveas(gcf,[outname,'.jpg'],'jpg');
% save([outname,'_obj.mat'],'CGobj');
% 
% % % removing non-finite values from C is same as removing zero-values from
% % % COV : run the lines below to verify this
% % C = NadjSPC./COV
% % [cx,cy] = find(~isfinite(C));
% % Cnf = [];
% % for i=1:length(cx), Cnf=[Cnf; C(cx(i),cy(i))]; end
% % [zx,zy] = find(COV==0);
% % COVz = [];
% % for i=1:length(zx), COVz = [COVz; COV(zx(i),zy(i))]; end
% 
% % ----------------------------------------------------

