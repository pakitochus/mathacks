function [performance,confMat] = structConfMat(perf)
% Returns a performance structure, with the values of Correct Rate, Error
% Rate, Sensitivity, Specificity, npv, ppv, Positive Likelihood and
% Negative Likelihood. It uses a 2x2 confusion matrix, so no higher
% dimension or non-binary classifiers are allowed. 

    confMat = sum(perf);
    confMat = reshape(confMat,[sqrt(numel(confMat)),sqrt(numel(confMat))]);
    performance = struct();
    N = sum(confMat(:)); 
    performance.CorrectRate = sum(diag(confMat))/N;
    performance.CRstd = std(sum(perf(:,[1,4]),2)./sum(perf,2));
    performance.ErrorRate = 1-performance.CorrectRate;
%     aux=diag(confMat)./sum(confMat,2);
    aux = perf(:,[1,end])./[sum(perf(:,[1,3]),2),sum(perf(:,[2,4]),2)];
    performance.Sensitivity = mean(aux(:,end)); % aux(end);
    performance.SensSTD = std(aux(:,end));
    performance.Specificity = mean(aux(:,1)); % aux(1);
    performance.SpecSTD = std(aux(:,1));
    aux=diag(confMat)./sum(confMat,1)';
    performance.ppv = aux(end);
    performance.npv = aux(1);
    performance.PositiveLikelihood = performance.Sensitivity/(1-performance.Specificity);
    performance.NegativeLikelihood = (1-performance.Sensitivity)/performance.Specificity; 
    performance.F1 = 2*performance.ppv*performance.Sensitivity/(performance.ppv+performance.Sensitivity);
    performance.precision = performance.Sensitivity/(performance.Sensitivity+1-performance.Specificity);
end

% Estructura en Matlab:
%          |Test Negative | Test Positive |
% -----------------------------------------
% Negative |     TN       |      FP       | 
% -----------------------------------------
% Positive |     FN       |      TP       |
% -----------------------------------------