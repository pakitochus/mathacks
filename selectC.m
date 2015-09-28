function [bestc,bestg] = selectC(ytrain,XTRAIN,kernel)
% Selects optimum C value in a linear kernel (using MATLAB) or another SVM kernel (using LIBSVM).
bestcv = 0;
XTRAIN = sparse(XTRAIN);
log2g=0;bestg=0;
for log2c = -10:5,
    
    if(kernel==0)
%         cmd=['-t ', num2str(kernel),' -v 5 -c ', num2str(2^log2c), ' -q'];
%         cv = svmtrain(ytrain, XTRAIN, cmd);
%         if (cv >= bestcv),
%             bestcv = cv; bestc = 2^log2c; 
%         end
                cmd = [' -v 5  -c ', num2str(2^log2c), ' -q'];
                cv = train(ytrain, XTRAIN, cmd);
                if (cv >= bestcv),
                    bestcv = cv; bestc = 2^log2c;
                end
    else
        for log2g = -22:2,
            cmd=['-t ', num2str(kernel),' -v 5 -c ', num2str(2^log2c), ' -g ', num2str(2^log2g), ' -q'];
            cv = svmtrain(ytrain, XTRAIN, cmd);
            if (cv >= bestcv),
                bestcv = cv; bestc = 2^log2c; bestg = 2^log2g;
            end
        end
    end
    fprintf('%g %g %g (best c=%g, g=%g, rate=%g)\n', log2c, log2g, cv, bestc, bestg, bestcv);
end
end
