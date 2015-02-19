function [IDX,Z]=rankfeaturesInterface(data, group, criterion)
% RANKFEATURESINTERFACE provides an interface to integrate filtering
% methods other than the five included in RankFeatures, such as FDR or 
% the p-value estimating method using ANOVA1, and expanding the t-test
% using the two-tailed, two-sample t-test. 
%

originals={'wilcoxon', 'entropy', 'roc', 'bhattacharyya'};
new={'ttest','fdr','pvalue','none'};
[N,M]=size(data);
if(any(strcmp(criterion, originals)))
    [IDX,Z]=rankfeatures(data', group, 'Criterion',criterion);
    Z(isnan(Z))=0;
elseif(any(strcmp(criterion,new)))
    switch(criterion)
        case 'ttest'
            [~,~,~,stats]=ttest2(data(~group,:),data(group,:));
            Z=stats.tstat';
            Z(isnan(Z))=0;
            [~, IDX]= sort(abs(Z),'descend');
        case 'fdr'
            Z = fdr(data, group)';
            Z(isnan(Z))=0;
            [~, IDX]= sort(Z,'descend');
        case 'none'
            IDX=1:M;
            %             IDX=size(data,2):-1:1;
        case 'pvalue'
            count=0;
            for comp=1:M
            	Z(comp) = anova1(data(:,comp),group,'off');
                if(mod(comp,10)==0)
                    fprintf(1, repmat('\b',1,count)); %delete line before
                    count=fprintf('%d%%',round(100*comp/M));
                end
            end
            Z(isnan(Z))=0;
            [~, IDX]= sort(abs(Z),'descend');
            fprintf('\n');
    end
else
    error('Valid criterion not introduced');
end
end

function [ fdr ] = fdr(data,labels)
%FDR compute point by point fisher discriminant ratio
% [ fdr ] = fdr(data,labels)
%  data: array of 2 dimentions, fdr is computed along the first dimention
%   labels: binary labels
% note: computations are tolerant to NANs

%compute FDR
m0=nanmean(data(find(labels),:),1);%means for sane group
v0=nanvar(double(data(find(labels),:)),0,1);%variance for sane group
m1=nanmean(data(find(~labels),:),1);%means for AD group
v1=nanvar(double(data(find(~labels),:)),0,1);%covariance for AD group
fdr=((m0-m1)).^2./(v1+v0);%compute Fisher Discriminant Ratio

end

