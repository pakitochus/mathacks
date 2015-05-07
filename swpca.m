function [Xhat, varargout] = swpca(X,cat,d,varargin)
% SWPCA performs the Statistical Significance Weighted Principal
% Component Analysis on a dataset, given a set of categorical variables
% (cat), the effect vector d of the variables (containing +1 if we want to
% enhance the effect of the variable and -1 if we want to reduce the effect
% of the variable) and k, the inverse of the significance threshold. It
% returns the enhanced dataset.
% USAGE:
%   Xhat = SWPCA(X,cat,d). computes the signal using k=0.05 and returns the
%   reconstructed signal. 
% 
%   Xhat = SWPCA(X,cat,d,'k',k) specifies the current k. 
% 
%   Xhat = SWPCA(X,cat,d,'training',trSet) gives the TRSET boolean vector
%   with the elements belonging to the training set. 
% 
%   [Xhat, S, Lambda, A] = SWPCA(...) returns the reconstructed signal
%   Xhat, the component scores S, the diagonal weighting matrix Lambda and
%   the pseudoinverse of the weightning matrix A. 

% Sets the defaults.
k = 0.05;
N = size(X,1);
trSet = true(N,1); % binary vector with the training set.
training = false;

if ~isempty(varargin)
    for i=1:length(varargin)
        parametro = varargin{i};
        if(ischar(parametro))
            switch(parametro)
                case 'k'
                    if (length(varargin)>i)&&(~isempty(varargin{i+1}))
                        k = varargin{i+1};
                    end
                case 'training'
                    if (length(varargin)>i)&&(~isempty(varargin{i+1}))
                        trSet = varargin{i+1};
                        training = true;
                    end
            end
        end
        
    end
end

% CHECK AND INIZIALIZE
if(size(cat,2)~=length(d))
    error('Number of columns in CAT must be equal to length of D');
end
meanTr = mean(X,2);
X = bsxfun(@minus,X,meanTr);
varTr = var(X,0,2);
X = bsxfun(@rdivide,X,varTr);
if(training)
    XNTEST = X(~trSet,:);
end
XNTRAIN = X(trSet,:);
clear X;

% PERFORM PCA OVER THE TRAINING SET
[W, ~]=princomp(XNTRAIN,'econ');
STr=XNTRAIN*W;

if(training)
    STe=XNTEST*W;
end

% PERFORM ANOVA ON THE COMPONENTS OF THE TRAINING SET
vars = cell(1,length(d));
for v=1:length(d)
    if(training)
        vars{v} = cat(trSet,v);
    else
        vars{v} = cat(:, v);
    end
end
C = size(STr,2);
p = zeros(C,length(d));
for comp=1:C
    p(comp,:)=anovan(STr(:,comp),vars,'model','linear','display','off');
end

% COMPUTE THE WEIGHTINGS
lambda=zeros(1,C);
for c = 1:C,
    lambda(c)=weighting(p(c,:)',d,k);
end

% RECONSTRUCT THE SIGNALS.
A = pinv(W);
Lambda = diag(lambda);
Xhat = STr*Lambda*A;
if(training)
    XTRhat = Xhat;
    XTEhat = STe*Lambda*A;
    Xhat = zeros(length(trSet),size(XTRhat,2));
    Xhat(trSet,:) = XTRhat;
    Xhat(~trSet,:) = XTEhat;
    clear XTRhat XTEhat;
end
Xhat = bsxfun(@times, Xhat, varTr);
Xhat = bsxfun(@plus, Xhat, meanTr);

% RETURN THE PARAMETERS
nout = max(nargout,1) - 1;
if nout>1
    if(training)
        S = zeros(N,L);
        S(training,:) = STr;
        S(~training,:) = STe;
        clear STr STe;
        varargout{1} = S;
    else
        varargout{1} = STr;
    end
end
if nout>2
    varargout{2} = Lambda;
end
if nout>3
    varargout{3} = A;
end
end


function y = weighting(p,d,k)
sumatoria = 1;
for jarl=1:length(d),
    if ndims(p)==1
        sumatoria = sumatoria + power(-1,1-d(jarl))*exp(-squeeze(abs(p(jarl)))/k);
    elseif ismatrix(p)
        sumatoria = sumatoria + power(-1,1-d(jarl))*exp(-squeeze(abs(p(jarl,:)))/k);
    elseif ndims(p)==3
        sumatoria = sumatoria + power(-1,1-d(jarl))*exp(-squeeze(abs(p(jarl,:,:)))/k);
    else
        error('P matrices with dimensions > 3 not supported');
    end
end
y = sumatoria/2; 
end
