function y = matrixToLatex(matrix,format)
% MATRIXTOLATEX(matrix) converts a 1D, 2D or 3D matrix into a latex table
% that can be imported into whichever latex document. 
matrix=squeeze(matrix);
tamano=size(matrix);
dims=ndims(matrix),
if(dims==1)
    for i=1:numel(matrix)
        fprintf(['$',format,'$ \\\\\n'], matrix(i));
    end
elseif(dims==2)
    for i=1:tamano(1)
        for j=1:tamano(2)
            if(j<tamano(2))
                fprintf(['$', format,'$ & '], matrix(i,j));
            else
                fprintf(['$',format,'$ '], matrix(i,j));
            end
        end
        fprintf('\\\\\n');
    end
elseif(dims==3)
    for i=1:tamano(1)
        for j=1:tamano(2)
            for k=1:tamano(3)
                if(j<tamano(3))
                    fprintf(['$', format,'$ & '], matrix(i,j,k));
                else
                    fprintf(['$',format,'$ '], matrix(i,j,k));
                end
            end
                fprintf('\\\\\n');
        end
            fprintf('\\hline\n');
    end
else
    error('matrix has more than 3 or less than 1 dimension');
end
    
