% Answering questions marked in yellow on the draft manuscript
clear

% -------------------------------------------------------------------------
% Q: Can we find genes with robust mRNA and protein data that show strongly negative correlations along the
% leaf gradient??? Most individual genes in protein/mRNA clusters show a general positive correlation, but
% a developmental shift in ratio. We don’t really have examples of inverse mRNA-protein abundances? Can
% we look for them?
% A: Since "robust" data is needed, let's look at the filtered dataset
% -------------------------------------------------------------------------
[N,T] = xlsread('Supplemental Table 1A-B-C-D-E-aug2013-v1.xls',5);
NSAF = N(:,4:7);
RPKM = N(:,8:11);
NadjSPC = N(:,13:16);
COV = N(:,17:20);
groupID = N(:,1);

% -- raw values --
% fprintf('\n\ngroupID\tNadjSPC_1\tNadjSPC_4\tNadjSPC_9\tNadjSPC_14\tCOV_1\tCOV_4\tCOV_9\tCOV_14\tCorr\tPvalue\n');
% for i=1:size(NadjSPC,1)
%     [r,p]=corr(NadjSPC(i,:)',COV(i,:)');
%     if (r<0)&&(p<0.1)
%         fprintf('%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n',groupID(i),NadjSPC(i,1),NadjSPC(i,2),NadjSPC(i,3),NadjSPC(i,4),COV(i,1),COV(i,2),COV(i,3),COV(i,4),r,p);
%         % plot(COV(i,:),NadjSPC(i,:),'b*'); lsline; pause
%     end
% end
fprintf('\ngroupID\tNSAF_1\tNSAF_4\tNSAF_9\tNSAF_14\tRPKM_1\tRPKM_4\tRPKM_9\tRPKM_14\tCorr\tPvalue\n');
for i=1:size(NSAF,1)
    [r,p]=corr(NSAF(i,:)',RPKM(i,:)');
    if (r<0)&&(p<0.1)
        fprintf('%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n',groupID(i),NSAF(i,1),NSAF(i,2),NSAF(i,3),NSAF(i,4),RPKM(i,1),RPKM(i,2),RPKM(i,3),RPKM(i,4),r,p);
    end
end

% -- log values with 1e-6 substituted in place of prot=0 --
% fprintf('\n\ngroupID\tNadjSPC_1\tNadjSPC_4\tNadjSPC_9\tNadjSPC_14\tCOV_1\tCOV_4\tCOV_9\tCOV_14\tCorr\tPvalue\n');
% for i=1:size(NadjSPC,1)
%     % ignore if any of the COV points are zero
%     if any(COV(i,:)==0), continue, end
%     % replace any prot=0 values by 1e-6
%     NadjSPC(i,find(NadjSPC(i,:)==0)) = 1e-6;
%     
%     [r,p]=corr(log10(NadjSPC(i,:)'),log10(COV(i,:)'));
%     if (r<0)&&(p<0.1)
%         fprintf('%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n',groupID(i),log10(NadjSPC(i,1)),log10(NadjSPC(i,2)),log10(NadjSPC(i,3)),log10(NadjSPC(i,4)),log10(COV(i,1)),log10(COV(i,2)),log10(COV(i,3)),log10(COV(i,4)),r,p);
%         % plot(COV(i,:),NadjSPC(i,:),'b*'); lsline; pause
%     end
% end
fprintf('\ngroupID\tNSAF_1\tNSAF_4\tNSAF_9\tNSAF_14\tRPKM_1\tRPKM_4\tRPKM_9\tRPKM_14\tCorr\tPvalue\n');
for i=1:size(NSAF,1)
    % ignore if any of the RPKM points are zero
    if any(RPKM(i,:)==0), continue, end
    % replace any prot=0 values by 1e-6
    NSAF(i,find(NSAF(i,:)==0)) = 1e-6;    

    [r,p]=corr(log10(NSAF(i,:)'),log10(RPKM(i,:)'));
    if (r<0)&&(p<0.1)
        fprintf('%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n',groupID(i),log10(NSAF(i,1)),log10(NSAF(i,2)),log10(NSAF(i,3)),log10(NSAF(i,4)),log10(RPKM(i,1)),log10(RPKM(i,2)),log10(RPKM(i,3)),log10(RPKM(i,4)),r,p);
        % plot(RPKM(i,:),NSAF(i,:),'b*'); lsline; pause
    end
end

% -------------------------------------------------------------------------


% % -------------------------------------------------------------------------
% % Q: The 1% of proteins (31) not detected at the mRNA level were.. and had all low abundance??
% % In the full matched dataset, find genes having mRNA=0 and examine their
% % protein abundance
% % -------------------------------------------------------------------------
% [N,T] = xlsread('forLalit-SupplTable1D-correlations.xlsx');
% NSAF = N(:,5:8);
% NadjSPC = N(:,9:12);
% COV = N(:,14:17);
% RPKM = N(:,18:21);
% groupID = N(:,1);
% 
% sec = {'sec1','sec4','sec9','sec14'};
% 
% % NadjSPC & COV
% [ix,iy] = find(COV==0);
% fprintf('GroupID\tSection\tNadjSPC\tCOV\n');
% for k=1:length(ix)
%     fprintf('%d\t%s\t%1.10f\t%1.10f\n',groupID(ix(k)),sec{iy(k)},NadjSPC(ix(k),iy(k)),COV(ix(k),iy(k)));
% end
% 
% % NSAF & RPKM
% [ix,iy] = find(RPKM==0);
% fprintf('GroupID\tSection\tNSAF\tRPKM\n');
% for k=1:length(ix)
%     fprintf('%d\t%s\t%1.10f\t%1.10f\n',groupID(ix(k)),sec{iy(k)},NSAF(ix(k),iy(k)),RPKM(ix(k),iy(k)));
% end
% % -------------------------------------------------------------------------


% % -------------------------------------------------------------------------
% % Q: Are there any genes that have high mRNA abundance but were not detected at the protein level??
% % Read the full matched dataset, find genes that had NSAF or NadjSPC=0 and
% % look at their mRNA abundances
% % -------------------------------------------------------------------------
% [N,T] = xlsread('forLalit-SupplTable1D-correlations.xlsx');
% NSAF = N(:,5:8);
% NadjSPC = N(:,9:12);
% COV = N(:,14:17);
% RPKM = N(:,18:21);
% groupID = N(:,1);
% 
% p = 0.9; % quantile
% 
% % NadjSPC & COV
% Y=[];
% [ix,iy]=find(NadjSPC==0);
% for k=1:length(ix)
%     if COV(ix(k),iy(k))>=quantile(COV(:),p)
%         Y = [Y; [groupID(ix(k)) ix(k) iy(k)]];
%     end
% end
% sec = {'sec1','sec4','sec9','sec14'};
% fprintf('GroupID\tSection\tNadjSPC\tCOV\n');
% for i=1:size(Y,1)
%     fprintf('%d\t%s\t%f\t%f\n',Y(i,1),sec{Y(i,3)},NadjSPC(Y(i,2),Y(i,3)),COV(Y(i,2),Y(i,3)));
% end
% 
% % NSAF & RPKM
% Y=[];
% [ix,iy]=find(NSAF==0);
% for k=1:length(ix)
%     if RPKM(ix(k),iy(k))>=quantile(RPKM(:),p)
%         Y = [Y; [groupID(ix(k)) ix(k) iy(k)]];
%     end
% end
% sec = {'sec1','sec4','sec9','sec14'};
% fprintf('GroupID\tSection\tNSAF\tRPKM\n');
% for i=1:size(Y,1)
%     fprintf('%d\t%s\t%f\t%f\n',Y(i,1),sec{Y(i,3)},NSAF(Y(i,2),Y(i,3)),RPKM(Y(i,2),Y(i,3)));
% end
% % -------------------------------------------------------------------------

