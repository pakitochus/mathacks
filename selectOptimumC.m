function [bestC,bestacc] = selectOptimumC(XTRAIN,ytrain)
    % Selection of optimum C by x-validation 
	% Requires: LIBSVM
    c=logspace(-2, 0, 13);
    bestacc=0;
    bestC=0;
    subcp = cvpartition(ytrain,'k',5);%Hacemos las partes
    for perc=1:numel(c),
        
    fprintf('.')
        C=c(perc);
        save('paramC.mat', 'C');
        performance = crossval(@featuresAsFeatures,XTRAIN,ytrain,'partition',subcp);
        perf = structConfMat(performance);
        acc(perc)=perf.CorrectRate;
        if(acc(perc)>bestacc+exp(perc/20)/150),
            bestC=c(perc);
            bestacc=acc(perc);
        end
    end
end

function testval = featuresAsFeatures(XTRAIN,ytrain,XTEST,ytest)

    if(exist('paramC.mat', 'file')==2)
        load('paramC.mat')
    else
        C=1;
    end
% options=optimset('MaxIter', 100000);
    % Con CLASSIFY linear, 
        if(size(XTRAIN,2)>0)
            try
                model = svmtrain(double(ytrain),XTRAIN,sprintf('-s 0 -t 0 -c %f -q',C));
            catch err
                error(err.message);
            end
            ypredicted=logical(svmpredict(double(ytest), XTEST, model,'-q'));
            testval = confusionmat(ytest,ypredicted, 'order', unique(ytrain));
        else
            return;
        end
end

