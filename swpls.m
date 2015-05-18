function [Xhat, varargout] = swpls(X,catVars,d,varargin)
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
if(size(catVars,2)~=length(d))
    error('Number of columns in CAT must be equal to length of D');
end

% PERFORM PLS ON THE TRAINING SET
if size(catVars,1)~=size(X,1),
    catVars = catVars';
end
meanTr = mean(X(trSet,:),1);
meanY = mean(catVars(trSet,:),1);
X0 = bsxfun(@minus, X, meanTr);
Y0 = bsxfun(@minus, catVars, meanY);
clear X; 

ncomp = 10; 
[W,~] = simpls(X0(trSet,:),Y0(trSet),ncomp);
STr = X0(trSet,:)*W; 
if(training)
    STe = X0(~trSet,:)*W; 
end
% meanTr = mean(X,2);
% X = bsxfun(@minus,X,meanTr);
% varTr = var(X,0,2);
% X = bsxfun(@rdivide,X,varTr);
% if(training)
%     XNTEST = X(~trSet,:);
% end
% XNTRAIN = X(trSet,:);
% clear X;

% % PERFORM PCA OVER THE TRAINING SET
% [W, ~]=princomp(XNTRAIN,'econ');
% STr=XNTRAIN*W;

% if(training)
%     STe=XNTEST*W;
% end

% A = stats.W;

% PERFORM ANOVA ON THE COMPONENTS OF THE TRAINING SET
vars = cell(1,length(d));
for v=1:length(d)
    if(training)
        vars{v} = catVars(trSet,v);
    else
        vars{v} = catVars(:, v);
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
% Xhat = bsxfun(@times, Xhat, varTr);
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
        sumatoria = sumatoria + d(jarl)*exp(-squeeze(abs(p(jarl)))/k);
    elseif ismatrix(p)
        sumatoria = sumatoria + d(jarl)*exp(-squeeze(abs(p(jarl,:)))/k);
    elseif ndims(p)==3
        sumatoria = sumatoria + d(jarl)*exp(-squeeze(abs(p(jarl,:,:)))/k);
    else
        error('P matrices with dimensions > 3 not supported');
    end
end
y = sumatoria/length(d);
end

%SIMPLS Basic SIMPLS.  Performs no error checking.
function [Xloadings,Yloadings,Xscores,Yscores,Weights] = simpls(X0,Y0,ncomp)

[n,dx] = size(X0);
dy = size(Y0,2);

% Preallocate outputs
outClass = superiorfloat(X0,Y0);
Xloadings = zeros(dx,ncomp,outClass);
Yloadings = zeros(dy,ncomp,outClass);
if nargout > 2
    Xscores = zeros(n,ncomp,outClass);
    Yscores = zeros(n,ncomp,outClass);
    if nargout > 4
        Weights = zeros(dx,ncomp,outClass);
    end
end

% An orthonormal basis for the span of the X loadings, to make the successive
% deflation X0'*Y0 simple - each new basis vector can be removed from Cov
% separately.
V = zeros(dx,ncomp);

Cov = X0'*Y0;
for i = 1:ncomp
    % Find unit length ti=X0*ri and ui=Y0*ci whose covariance, ri'*X0'*Y0*ci, is
    % jointly maximized, subject to ti'*tj=0 for j=1:(i-1).
    [ri,si,ci] = svd(Cov,'econ'); ri = ri(:,1); ci = ci(:,1); si = si(1);
    ti = X0*ri;
    normti = norm(ti); ti = ti ./ normti; % ti'*ti == 1
    Xloadings(:,i) = X0'*ti;
    
    qi = si*ci/normti; % = Y0'*ti
    Yloadings(:,i) = qi;
    
    if nargout > 2
        Xscores(:,i) = ti;
        Yscores(:,i) = Y0*qi; % = Y0*(Y0'*ti), and proportional to Y0*ci
        if nargout > 4
            Weights(:,i) = ri ./ normti; % rescaled to make ri'*X0'*X0*ri == ti'*ti == 1
        end
    end
    
    % Update the orthonormal basis with modified Gram Schmidt (more stable),
    % repeated twice (ditto).
    vi = Xloadings(:,i);
    for repeat = 1:2
        for j = 1:i-1
            vj = V(:,j);
            vi = vi - (vj'*vi)*vj;
        end
    end
    vi = vi ./ norm(vi);
    V(:,i) = vi;
    
    % Deflate Cov, i.e. project onto the ortho-complement of the X loadings.
    % First remove projections along the current basis vector, then remove any
    % component along previous basis vectors that's crept in as noise from
    % previous deflations.
    Cov = Cov - vi*(vi'*Cov);
    Vi = V(:,1:i);
    Cov = Cov - Vi*(Vi'*Cov);
end

if nargout > 2
    % By convention, orthogonalize the Y scores w.r.t. the preceding Xscores,
    % i.e. XSCORES'*YSCORES will be lower triangular.  This gives, in effect, only
    % the "new" contribution to the Y scores for each PLS component.  It is also
    % consistent with the PLS-1/PLS-2 algorithms, where the Y scores are computed
    % as linear combinations of a successively-deflated Y0.  Use modified
    % Gram-Schmidt, repeated twice.
    for i = 1:ncomp
        ui = Yscores(:,i);
        for repeat = 1:2
            for j = 1:i-1
                tj = Xscores(:,j);
                ui = ui - (tj'*ui)*tj;
            end
        end
        Yscores(:,i) = ui;
    end
end
end
