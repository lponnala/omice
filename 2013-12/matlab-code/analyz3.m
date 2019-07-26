% // By Function //
clear;

% ----------------------------------------------------------------------
% "Function X Section" average abundances
% ----------------------------------------------------------------------
[N,T,R] = xlsread('forLalit-SupplTable1E-by-plastidorfunctionv2.xlsx');
NSAF = N(:,10:13);
RPKM = N(:,14:17);
NadjSPC = N(:,19:22);
COV = N(:,23:26);

func = cell(size(R,1)-1,1);
for i=2:size(R,1)
    elem = R{i,3};
    if isnumeric(elem)
        if ~isnan(elem)
            func{i-1} = num2str(elem);
        else
            func{i-1} = '';
        end
    elseif ischar(elem)
        func{i-1} = elem;
    end
end

sec = {'sec1','sec4','sec9','sec14'};

fid = unique(func);
fprintf('\nUsing all proteins\n');
fprintf('Section\t#Proteins\tavg_NadjSPC\tavg_COV\tavg_NSAF\tavg_RPKM\n')
for i=1:length(fid)
    if strcmp(fid{i},''), continue, end
    fprintf('Function: %s\n',fid{i});
    for j=1:length(sec)
        ind = find(strcmp(func,fid{i}));
        fprintf('%s\t%d\t%f\t%f\t%f\t%f\n',sec{j},length(ind),mean(NadjSPC(ind,j)),mean(COV(ind,j)),mean(NSAF(ind,j)),mean(RPKM(ind,j)));
    end
end
fprintf('\nExcluding proteins that have NadjSPC=0\n');
fprintf('Section\t#Proteins\tavg_NadjSPC\tavg_COV\tavg_NSAF\tavg_RPKM\n')
for i=1:length(fid)
    if strcmp(fid{i},''), continue, end
    fprintf('Function: %s\n',fid{i});
    for j=1:length(sec)
        ind = find(strcmp(func,fid{i}) & NadjSPC(:,j)>0);
        fprintf('%s\t%d\t%f\t%f\t%f\t%f\n',sec{j},length(ind),mean(NadjSPC(ind,j)),mean(COV(ind,j)),mean(NSAF(ind,j)),mean(RPKM(ind,j)));
    end
end
% ----------------------------------------------------------------------


% % ----------------------------------------------------------------------
% % Function-wise correlation within each section, leave out zero values,
% % with log transformation
% % ----------------------------------------------------------------------
% [N,T,R] = xlsread('forLalit-SupplTable1E-by-plastidorfunctionv2.xlsx');
% NSAF = N(:,10:13);
% RPKM = N(:,14:17);
% NadjSPC = N(:,19:22);
% COV = N(:,23:26);
% 
% func = cell(size(R,1)-1,1);
% for i=2:size(R,1)
%     elem = R{i,3};
%     if isnumeric(elem)
%         if ~isnan(elem)
%             func{i-1} = num2str(elem);
%         else
%             func{i-1} = '';
%         end
%     elseif ischar(elem)
%         func{i-1} = elem;
%     end
% end
% 
% sec = {'sec1','sec4','sec9','sec14'};
% 
% fid = unique(func);
% fprintf('\nNadjSPC & COV\n')
% fprintf('Section\t#Proteins\tRp\tp.value\tRs\tp.value\n')
% for i=1:length(fid)
%     if strcmp(fid{i},''), continue, end
%     fprintf('Function: %s\n',fid{i});
%     for j=1:length(sec)
%         ind = find(strcmp(func,fid{i}) & NadjSPC(:,j)>0);
%         [Rp,Pp] = corr(log10(NadjSPC(ind,j)),log10(COV(ind,j)),'type','Pearson');
%         [Rs,Ps] = corr(log10(NadjSPC(ind,j)),log10(COV(ind,j)),'type','Spearman');
%         fprintf('%s\t%d\t%f\t%f\t%f\t%f\n',sec{j},length(ind),Rp,Pp,Rs,Ps);
%     end
% end
% fprintf('\nNSAF & RPKM\n')
% fprintf('Section\t#Proteins\tRp\tp.value\tRs\tp.value\n')
% for i=1:length(fid)
%     if strcmp(fid{i},''), continue, end
%     fprintf('Function: %s\n',fid{i});
%     for j=1:length(sec)
%         ind = find(strcmp(func,fid{i}) & NadjSPC(:,j)>0);
%         [Rp,Pp] = corr(log10(NSAF(ind,j)),log10(RPKM(ind,j)),'type','Pearson');
%         [Rs,Ps] = corr(log10(NSAF(ind,j)),log10(RPKM(ind,j)),'type','Spearman');
%         fprintf('%s\t%d\t%f\t%f\t%f\t%f\n',sec{j},length(ind),Rp,Pp,Rs,Ps);
%     end
% end
% % ----------------------------------------------------------------------


% % ----------------------------------------------------------------------
% % Function-wise correlation within each section, leave out zero values,
% % no log transformation
% % ----------------------------------------------------------------------
% [N,T,R] = xlsread('forLalit-SupplTable1E-by-plastidorfunctionv2.xlsx');
% NSAF = N(:,10:13);
% RPKM = N(:,14:17);
% NadjSPC = N(:,19:22);
% COV = N(:,23:26);
% 
% func = cell(size(R,1)-1,1);
% for i=2:size(R,1)
%     elem = R{i,3};
%     if isnumeric(elem)
%         if ~isnan(elem)
%             func{i-1} = num2str(elem);
%         else
%             func{i-1} = '';
%         end
%     elseif ischar(elem)
%         func{i-1} = elem;
%     end
% end
% 
% sec = {'sec1','sec4','sec9','sec14'};
% 
% fid = unique(func);
% fprintf('\nNadjSPC & COV\n')
% fprintf('Section\t#Proteins\tRp\tp.value\tRs\tp.value\n')
% for i=1:length(fid)
%     if strcmp(fid{i},''), continue, end
%     fprintf('Function: %s\n',fid{i});
%     for j=1:length(sec)
%         ind = find(strcmp(func,fid{i}) & NadjSPC(:,j)>0);
%         [Rp,Pp] = corr(NadjSPC(ind,j),COV(ind,j),'type','Pearson');
%         [Rs,Ps] = corr(NadjSPC(ind,j),COV(ind,j),'type','Spearman');
%         fprintf('%s\t%d\t%f\t%f\t%f\t%f\n',sec{j},length(ind),Rp,Pp,Rs,Ps);
%     end
% end
% fprintf('\nNSAF & RPKM\n')
% fprintf('Section\t#Proteins\tRp\tp.value\tRs\tp.value\n')
% for i=1:length(fid)
%     if strcmp(fid{i},''), continue, end
%     fprintf('Function: %s\n',fid{i});
%     for j=1:length(sec)
%         ind = find(strcmp(func,fid{i}) & NadjSPC(:,j)>0);
%         [Rp,Pp] = corr(NSAF(ind,j),RPKM(ind,j),'type','Pearson');
%         [Rs,Ps] = corr(NSAF(ind,j),RPKM(ind,j),'type','Spearman');
%         fprintf('%s\t%d\t%f\t%f\t%f\t%f\n',sec{j},length(ind),Rp,Pp,Rs,Ps);
%     end
% end
% % ----------------------------------------------------------------------


% % ----------------------------------------------------------------------
% % Function-wise correlation within each section, no cutoff
% % ----------------------------------------------------------------------
% [N,T,R] = xlsread('forLalit-SupplTable1E-by-plastidorfunctionv2.xlsx');
% NSAF = N(:,10:13);
% RPKM = N(:,14:17);
% NadjSPC = N(:,19:22);
% COV = N(:,23:26);
% 
% func = cell(size(R,1)-1,1);
% for i=2:size(R,1)
%     elem = R{i,3};
%     if isnumeric(elem)
%         if ~isnan(elem)
%             func{i-1} = num2str(elem);
%         else
%             func{i-1} = '';
%         end
%     elseif ischar(elem)
%         func{i-1} = elem;
%     end
% end
% 
% sec = {'sec1','sec4','sec9','sec14'};
% 
% fid = unique(func);
% fprintf('\nNadjSPC & COV\n')
% fprintf('Section\t#Proteins\tRp\tp.value\tRs\tp.value\n')
% for i=1:length(fid)
%     if strcmp(fid{i},''), continue, end
%     fprintf('Function: %s\n',fid{i});
%     ind = find(strcmp(func,fid{i}));
%     for j=1:length(sec)
%         [Rp,Pp] = corr(NadjSPC(ind,j),COV(ind,j),'type','Pearson');
%         [Rs,Ps] = corr(NadjSPC(ind,j),COV(ind,j),'type','Spearman');
%         fprintf('%s\t%d\t%f\t%f\t%f\t%f\n',sec{j},length(ind),Rp,Pp,Rs,Ps);
%     end
% end
% fprintf('\nNSAF & RPKM\n')
% fprintf('Section\t#Proteins\tRp\tp.value\tRs\tp.value\n')
% for i=1:length(fid)
%     if strcmp(fid{i},''), continue, end
%     fprintf('Function: %s\n',fid{i});
%     ind = find(strcmp(func,fid{i}));
%     for j=1:length(sec)
%         [Rp,Pp] = corr(NSAF(ind,j),RPKM(ind,j),'type','Pearson');
%         [Rs,Ps] = corr(NSAF(ind,j),RPKM(ind,j),'type','Spearman');
%         fprintf('%s\t%d\t%f\t%f\t%f\t%f\n',sec{j},length(ind),Rp,Pp,Rs,Ps);
%     end
% end
% % ----------------------------------------------------------------------


% % ----------------------------------------------------------------------
% % -- in response to email from Klaas on 9/17/2013 --
% % ----------------------------------------------------------------------
% [N,T] = xlsread('lalit-forfigure3.xlsx');
% NSAF = N(:,10:13);
% RPKM = N(:,14:17);
% NadjSPC = N(:,20:23);
% COV = N(:,24:27);
% 
% func = T(2:end,3);
% sec = {'sec1','sec4','sec9','sec14'};
% cutoff = 0.00005;
% 
% fid = unique(func);
% fprintf('NadjSPC\tCOV\tNSAF\tRPKM\n');
% for i=1:length(fid)
%     fprintf('\n%s\n',fid{i});
%     for j=1:4
%         ind = find(strcmp(fid{i},func) & NadjSPC(:,j)>=cutoff);
%         fprintf('%s\t%d\t%f\t%f\t%f\t%f\n',sec{j},length(ind),mean(NadjSPC(ind,j)),mean(COV(ind,j)),mean(NSAF(ind,j)),mean(RPKM(ind,j)));
%     end
% end
% % ----------------------------------------------------------------------


% % ----------------------------------------------------------------------
% [N,T,R] = xlsread('forLalit-SupplTable1E-by-plastidorfunctionv2.xlsx');
% NSAF = N(:,10:13);
% RPKM = N(:,14:17);
% NadjSPC = N(:,19:22);
% COV = N(:,23:26);
% 
% func = cell(size(R,1)-1,1);
% for i=2:size(R,1)
%     elem = R{i,3};
%     if isnumeric(elem)
%         if ~isnan(elem)
%             func{i-1} = num2str(elem);
%         else
%             func{i-1} = '';
%         end
%     elseif ischar(elem)
%         func{i-1} = elem;
%     end
% end
% 
% fid = unique(func);
% fprintf('\nlog10(NadjSPC) & log10(COV)\n')
% fprintf('Function\t#Proteins\tRp\tp.value\tRs\tp.value\n')
% for i=1:length(fid)
%     if strcmp(fid{i},''), continue, end
%     ind = find(strcmp(func,fid{i}));
%     [Rp,Pp] = corr(log10(mean(NadjSPC(ind,:),2)),log10(mean(COV(ind,:),2)),'type','Pearson');
%     [Rs,Ps] = corr(log10(mean(NadjSPC(ind,:),2)),log10(mean(COV(ind,:),2)),'type','Spearman');
%     fprintf('%s\t%d\t%f\t%f\t%f\t%f\n',fid{i},length(ind),Rp,Pp,Rs,Ps);
% end
% fprintf('\nlog10(NSAF) & log10(RPKM)\n')
% fprintf('Function\t#Proteins\tRp\tp.value\tRs\tp.value\n')
% for i=1:length(fid)
%     if strcmp(fid{i},''), continue, end
%     ind = find(strcmp(func,fid{i}));
%     [Rp,Pp] = corr(log10(mean(NSAF(ind,:),2)),log10(mean(RPKM(ind,:),2)),'type','Pearson');
%     [Rs,Ps] = corr(log10(mean(NSAF(ind,:),2)),log10(mean(RPKM(ind,:),2)),'type','Spearman');
%     fprintf('%s\t%d\t%f\t%f\t%f\t%f\n',fid{i},length(ind),Rp,Pp,Rs,Ps);
% end
% % ----------------------------------------------------------------------


% % ---- OLD CODE ----
% [N,T,R] = xlsread('forLalit-SupplTable1E-by-simpleFunction.xlsx');
% NSAF = N(:,8:11);
% RPKM = N(:,12:15);
% NadjSPC = N(:,17:20);
% COV = N(:,21:24);
% 
% func = cell(size(R,1)-1,1);
% for i=2:size(R,1)
%     elem = R{i,3};
%     if isnumeric(elem)
%         if ~isnan(elem)
%             func{i-1} = num2str(elem);
%         else
%             func{i-1} = '';
%         end
%     elseif ischar(elem)
%         func{i-1} = elem;
%     end
% end
% 
% fid = unique(func);
% fprintf('\nNadjSPC & COV\n')
% for i=1:length(fid)
%     if strcmp(fid{i},''), continue, end
%     ind = find(strcmp(func,fid{i}));
%     [Rp,Pp] = corr(log10(mean(NadjSPC(ind,:),2)),log10(mean(COV(ind,:),2)),'type','Pearson');
%     [Rs,Ps] = corr(log10(mean(NadjSPC(ind,:),2)),log10(mean(COV(ind,:),2)),'type','Spearman');
%     fprintf('%s\t%d\t%f\t%f\t\t%f\t%f\n',fid{i},length(ind),Rp,Pp,Rs,Ps);
% end
% fprintf('\nNSAF & RPKM\n')
% for i=1:length(fid)
%     if strcmp(fid{i},''), continue, end
%     ind = find(strcmp(func,fid{i}));
%     [Rp,Pp] = corr(log10(mean(NSAF(ind,:),2)),log10(mean(RPKM(ind,:),2)),'type','Pearson');
%     [Rs,Ps] = corr(log10(mean(NSAF(ind,:),2)),log10(mean(RPKM(ind,:),2)),'type','Spearman');
%     fprintf('%s\t%d\t%f\t%f\t\t%f\t%f\n',fid{i},length(ind),Rp,Pp,Rs,Ps);
% end
