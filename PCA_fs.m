function [fs Rs]=PCA_fs(data,labels,npc,vethr)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Seleccion de características con PCA
% 
% A. Ortiz, 2014
% npc -> Número de componentes principales
% vethr -> umbral en la varianza acumulada (opcional)
% fs -> caracteristicas ordenadas por relevancia (primero la más relevante)
% Rs -> peso de cada caracteristica ordenadas por relevancia
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin <3
    error('Sintax: [fs]=PCA_fs(data,labels,npc,[vethr])');
end


%npc=2;
%nfeat=2;
%vethr=0.8;

lu=unique(labels);
[coef0 score0 lat0]=princomp(data(labels==lu(1),:)','econ');
[coef1 score1 lat1]=princomp(data(labels==lu(2),:)','econ');
ve0=cumsum(lat0)./sum(lat0);
ve1=cumsum(lat1)./sum(lat1);

if isempty(npc)
       npc0=find(ve0>vethr);
       npc1=find(ve1>vethr);
       npc=max([npc0(1) npc1(1)]);
       fprintf('npc no especificado... Usando %d PCs cubriendo el %d%% de la varianza...\n',npc,vethr*100);  
end


V0=[score0(:,1:npc)];
V1=[score1(:,1:npc)];

R=abs(max(abs(V0),[],2)-max(abs(V1),[],2));

[Rs fs]=sort(R,'descend');
%XX=[data(:,fs(1:nfeat))];
