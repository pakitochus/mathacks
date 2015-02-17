function [map,azim,elev] = mapBrain(I, varargin)
%MAPBRAIN performs a proyection of a 3D brain (from a struct I that
%contains an field img) into a two dimensional plane, using the spherical
%coordinates and different approaches.
%
%   [MAPA, COORD] = MAPBRAIN(I) takes the struct that contains the 3D image
%   in the field img, and projects it into two dimensions. The default
%   behaviour is a projection of the surface of the image I (using the
%   distance from the center to the last voxel with an intensity greater
%   than 0), just as in [1,2].
%
%   MAPA contains the 2D projected image using different approaches
%   (default, the surface of the image I).
%
%   COORD contains the 3D vectors that are used to project the image.
%
%   MAPBRAIN(..., 'resolution', R) sets the angle resolution to be used in
%   the projection. Default is 1 degree.
%
%   MAPBRAIN(..., 'deformation', D) sets the rate of unequally distributed
%   mapping vectors, to be used when the surface to be mapped is not
%   spherical but ellipsoid.
%
%   MAPBRAIN(..., 'threshold', TH) sets the threshold for the projections
%   needing it (surface, thickness, numfold, mediahigh).
%
%   MAPBRAIN(..., 'nlayers', N) sets up the layered projection, resulting
%   in a set of N projections, each one based on one of N divisions of the
%   r vector, depicting from the most internal to the most external regions
%
%   in MAPBRAIN(..., type), type is a string that specifies which operation
%   is applied to the mapping vectors to obtain the final projection. The
%   following projection types are defined in this work:
%       'sum' calculates the sum of all values in the mapping vector.
%       This is the default behaviour.
%
%       'surface' calculates the external surface of the tissue, using the
%       distance from the center of the image to the most external voxel
%       with an intensity greater than TH.
%
%       'thickness' measures the thickness of the tissue by calculating the
%       distance from the first voxel greater than a threshold TH to the
%       last voxel that acomplishes this criterion.
%
%       'average' performs the projection by averaging al the
%       values of intensity in the mapping vector.
%
%       'mediahigh' performs the projection by averaging al the
%       values of intensity greater than a threshold in the mapping vector.
%
%       'var' performs the projection by obtaining the variance
%       of the values of the intensity in the mapping vector.
%
%       'kurtosis' performs the projection by obtaining the kurtosis
%       of the values of the intensity in the mapping vector.
%
%       'skewness' performs the projection by obtaining the skewness
%       of the values of the intensity in the mapping vector.
%
%       'entropy' performs the projection by obtaining an estimate of the
%       entropy of the mapping vector.
%
%       'numfold' performs the projection by estimating the
%       number of cortex folds that the vector R crosses.
%
% created with MATLAB ver.: 7.14.0.739 (R2012a) on Linux 3.11.0-15-generic
% #25-Ubuntu SMP Thu Jan 30 17:25:07 UTC 2014 i686
%
% created by: fjesusmartinez@ugr.es
% UPDATED: Jan 9th 2015
%
% REFs:
% [1] - F.J. Martinez-Murcia et al. Projecting MRI Brain images for the
%       detection of Alzheimer's Disease. 2014.
% [2] - F.J. Martinez-Murcia et al. A Statistical Projection of MRI Brain
%       Images Approach for the Detection of Alzheimer's Disease. Journal
%       of Current Alzheimer's Disease N(V). 2015.

%Default Parameters:
res =1;                 % Angular Resolution
deformation=0;          % Deformation ratio
umbral=0;               % Threshold value
mostrarImagen = false;  % Checks if image is to be shown.
nlayers=1;              % The whole radius is considered. 

% Switches for different projections
projection='sum';

%Leemos varargin
if ~isempty(varargin)
    for i=1:length(varargin)
        parametro = varargin{i};
        if(ischar(parametro))
            switch(parametro)
                case 'resolution'
                    if (length(varargin)>i)&&(~isempty(varargin{i+1}))
                        res = varargin{i+1};
                    end
                case 'deformation'
                    if (length(varargin)>i)&&(~isempty(varargin{i+1}))
                        deformation = varargin{i+1};
                    end
                case 'threshold'
                    if (length(varargin)>i)&&(~isempty(varargin{i+1}))
                        th = varargin{i+1};
                        umbral = max(I.img(:))*th;
                    end
                case {'sum','thickness','surface','average','numfold',...
                        'mediahigh','var','skewness','entropy','kurtosis'}
                    projection = parametro;
                case 'nlayers'
                    if (length(varargin)>i)&&(~isempty(varargin{i+1}))
                        nlayers= varargin{i+1};
                    end
                case 'verbose'
                    if (length(varargin)>i)&&(~isempty(varargin{i+1}))
                        mostrarImagen = true;
                    end
                otherwise
                    error(['The option ',parametro,' is not recognized.']);
            end
        elseif(~ischar(parametro))
            if(~ischar(varargin{i-1}))
                error(['Not recognized option number ',num2str(i)]);
            end
        end
    end
end

% I.img = permute(I.img, [1 3 2]);
tam=size(I.img);
puntoMedio = ceil(size(I.img)/2);

% We create a set of vectors to map the coordinates.
spaceVector=1-deformation*cos(deg2rad(-3*180:res*2:180));
azim = deg2rad(cumsum(spaceVector)*res-270);
elev = deg2rad(90:-res:-90);
lon = length(0:max(puntoMedio));
tamArr=repmat(tam',1,lon);
[X,Y,Z]=meshgrid(azim,elev,0:max(puntoMedio));

[x,y,z] = sph2cart(X,Y,Z);
X=round(x+puntoMedio(1));
Y=round(y+puntoMedio(2));
Z=round(z+puntoMedio(3));

X(X>tamArr(1))=tamArr(1);X(X<1)=1;
Y(Y>tamArr(2))=tamArr(2);Y(Y<1)=1;
Z(Z>tamArr(3))=tamArr(3);Z(Z<1)=1;
coord=permute(sub2ind(tam,X,Y,Z), [2 1 3]);

% And map the selected vectors.
map = zeros(nlayers,size(coord,1), size(coord,2));
for nl=1:nlayers
    intvl=floor(lon/nlayers);
    for i=1:size(coord,1)
        for j=1:size(coord,2)
            radius=squeeze(I.img(coord(i,j,(1+(nl-1)*intvl):(nl*intvl))));
            switch(projection)
                case 'sum'
                    map(nl,i,j)=sum(radius(~isnan(radius)));
                case 'average'
                    map(nl,i,j)=mean(radius(~isnan(radius)));
                case 'var'
                    map(nl,i,j)=var(radius(~isnan(radius)));
                case 'skewness'
                    map(nl,i,j)=skewness(radius(~isnan(radius)));
                case 'kurtosis'
                    map(nl,i,j)=kurtosis(radius(~isnan(radius)));
                case 'entropy'
                    map(nl,i,j)=sum(radius((~isnan(radius)&(radius>0))).*log(radius((~isnan(radius)&(radius>0)))));
                case 'mediahigh'
                    map(nl,i,j)=mean(radius((radius>umbral)&(~isnan(radius))));
                case 'thickness'
                    radius=radius(~isnan(radius));
                    map(nl,i,j)=sum(radius>umbral);
                case 'numfold'
                    indices=find(radius>umbral);
                    if(numel(indices)>1)
                        t1=radius(indices(1):indices(end)-1)>umbral;
                        t2=radius(indices(1)+1:indices(end))>umbral;
                        map(nl,i,j)=sum(xor(t1,t2));
                    else
                        map(nl,i,j)=0;
                    end
                case 'surface'
                    if(isempty(find(radius>umbral, 1, 'last' )))
                        map(nl,i,j)=0;
                    else
                        map(nl,i,j)= find(radius>umbral, 1, 'last' );
                    end
            end
            %Check for changes
            %Returns the distance
        end
    end
end

if(mostrarImagen)
    set(0, 'defaultTextInterpreter', 'latex'); 
    for i=1:nlayers
        figure(i);
        imagesc(azim,elev,rot90(squeeze(map(i,:,:))));colormap('bone');
        xlabel('Azimuth $\Theta$'); 
        hl=ylabel('Elevation $\varphi$');%set(hl,'Interpreter','latex');
        axis image;
    end
end
map=squeeze(map);
end


