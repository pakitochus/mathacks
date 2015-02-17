function [bestc,bestg] = selectC(ytrain,XTRAIN,kernel)
bestcv = 0;
for log2c = -15:5,
    
    if(kernel==0)
        cmd = [' -v 5  -c ', num2str(2^log2c), ' -q'];
        cv = train(ytrain, sparse(XTRAIN), cmd);
        if (cv >= bestcv),
            bestcv = cv; bestc = 2^log2c;
        end
        log2g=NaN;
        bestg=NaN;
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
