function [performance,confMat] = structConfMatErr(perf)
% Returns a performance structure, with the values of Correct Rate, Error
% Rate, Sensitivity, Specificity, npv, ppv, Positive Likelihood and
% Negative Likelihood. It uses a 2x2 confusion matrix, so no higher
% dimension or non-binary classifiers are allowed. 
% Within each performance result, a second column with the standard error
% is added. 

    K = size(perf,1); % Number of folds
    N = sum(perf(:)); % Total number of individuals 
    Ncol = sum(perf,2); % Number of individuals per fold. 
    performance = struct();
    
    data = (perf(:,1)+perf(:,end))./Ncol; 
    performance.CorrectRate = [mean(data), std(data)/K];
    performance.ErrorRate = [mean(1-data), std(1-data)/K];
    
    data = perf(:,end)./sum(perf(:,2:2:end),2);
    performance.Sensitivity = [mean(data), std(data)/K];
    
    data = perf(:,1)./sum(perf(:,1:2:end),2);
    performance.Specificity = [mean(data), std(data)/K];
    
    data = perf(:,end)./sum(perf(:,end-1:end),2);
    performance.ppv = [mean(data), std(data)/K];
    
    data = perf(:,1)./sum(perf(:,1:2),2);
    performance.npv = [mean(data), std(data)/K];
    
    data = performance.Sensitivity(1)/(1-performance.Specificity(1));
    err = performance.Sensitivity(2)*(1/(1-performance.Specificity(1)))^2 + performance.Specificity(2)*(performance.Sensitivity(1)/(1-performance.Specificity(1))^2)^2;
    performance.PositiveLikelihood = [data, sqrt(err)/K];
    
    data=(1-performance.Sensitivity(1))/performance.Specificity(1); 
    err = performance.Sensitivity(2)*(1/(performance.Specificity(1)))^2 + performance.Specificity(2)*((1-performance.Sensitivity(1))/(performance.Specificity(1))^2)^2;
    performance.NegativeLikelihood = [data, sqrt(err)/K];
    
    performance.F1 = 2*performance.ppv(1)*performance.Sensitivity(1)/(performance.ppv(1)+performance.Sensitivity(1));
end

% Structure in Matlab:
%          |Test Negative | Test Positive |
% -----------------------------------------
% Negative |     TN       |      FP       | 
% -----------------------------------------
% Positive |     FN       |      TP       |
% -----------------------------------------