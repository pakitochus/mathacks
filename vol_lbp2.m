function [declbpfeat binlbpfeat]=vol_lbp2(I,P,R)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D LBP Implementation (Volumetric Local Binary Pattern)
% Spiral sampling is used, starting and ending at cylinder centre.
% generateRadialFilter from efficientLBP toolbox is required.
%
% Andrés Ortiz, January 2015.
% 
% 
% Ex: lbp=VOL_LBP(image,sampling_points,radius_from_center)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Debug
%I(:,:,1)=[100 142 100;131 130 125; 100 118 100];
%I(:,:,2)=[100 122 100;135 123 130; 100 120 100];
%I(:,:,3)=[100 134 100;123 129 129; 100 119 100];
% P=4,R=1 -> It should give 1  1  1  1  0  1  0  1  0  1  1  1  0  1 = 11951d

filtR=generateRadialFilterLBP(P, R);
[cx cy cz]=cube_center(I);
thrlevel=I(cx,cy,cz);
a=zeros(size(I,1),size(I,2));
ofs=size(a,1)-size(filtR,1);
filtR2=zeros(size(I,1),size(I,2),size(filtR,3));
if ofs>0
    a=zeros(size(I,1),size(I,2));
    a(:,1:ofs/2)=NaN;
    a(:,end-ofs/2+1)=NaN;
    a(1:ofs/2,:)=NaN;
    a(end-ofs/2+1,:)=NaN;
    b=a;
    for i=1:size(filtR,3)
        a=b;
        a(a==0)=filtR(:,:,i);
        filtR2(:,:,i)=a;
    end;
    filtR=filtR2;
    filtR(isnan(filtR))=0;
end;


binstr=[];
Ibin=zeros(size(I,1),size(I,2),size(I,2));
for k=1:size(I,3)   % k es la capa
    for i=1:P % i indica el vecino
        Ithr=(I(:,:,k)-thrlevel).*filtR(:,:,i);
        if (k==1)&&(i==1) % Primera capa y centro -> primer vecino
            % [cx cy]  ->  Valor central
            [x1 y1]=find(filtR(:,:,i)>0); %Siguiente valor;
            binstr=[binstr (Ithr(cx,cy).*filtR(cx,cy,i))>=0 (Ithr(x1,y1).*filtR(x1,y1,i))>=0]; 
        elseif k==size(I,3)&&(i==P)   % Ultima capa y último vecino  -> centro
            % [cx cy]  ->  Valor central
            if sum(sum(filtR(:,:,i)~=0)) > 2   % Valor Interpolado
                Ithr(cx,cy)=-Ithr(cx,cy);
                binstr=[binstr sum(Ithr(:))>=0];
            else    
                [x1 y1]=find(filtR(:,:,i)>0);  % Valor penúltimo
                binstr=[binstr (Ithr(x1,y1).*filtR(x1,y1,i))>=0 (Ithr(cx,cy).*filtR(cx,cy,i))>=0];
            end;
        else   
            % Para pixeles intermedios se multiplica por todo el filtro
            % para considerar valores interpolados
            if sum(sum(filtR(:,:,i)~=0)) > 2   % Valor Interpolado
                Ithr(cx,cy)=-Ithr(cx,cy);
                binstr=[binstr sum(Ithr(:))>=0];
            else
                [x1 y1]=find(filtR(:,:,i)>0); % Valor no central
                binstr=[binstr (Ithr(x1,y1)*filtR(x1,y1,i))>=0];
            end;
        end;
    end;
end;

A=(fliplr(binstr));
declbpfeat=sum(A.*2.^(numel(A)-1:-1:0));  % Conversion binario -> decimal
binlbpfeat=fliplr(binstr);



