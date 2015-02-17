function img = slicesDisplay(img4D,slices)
%   SLICESDISPLAY takes a 4-D image (with the first singleton dimension, 
%   such as the ones needed for the use within MONTAGE) and returns an 
%   image showing the selected slices side by side. 
verCorr=false;
numfilas=1;
if numel(slices)>9
    numfilas=ceil(numel(slices)/9);
    numcol=9;
else
    numcol=numel(slices);
end
tamx=size(img4D,2);
tamy=size(img4D,3);
img=zeros(tamx*numfilas,tamy*(numcol-1));
if(size(img4D,1)==1)
    if(size(squeeze(img4D(1,:,:,1)),1)<size(squeeze(img4D(1,:,:,1)),2)),
        verCorr=true;
    end
    for i=1:numel(slices),
        if(verCorr),
            img(((ceil(i/9)-1)*tamy+1):((ceil(i/9)-1)*tamy+tamy),(mod(i-1,9)*tamx+1):(mod(i-1,9)*tamx+tamx))= imrotate(squeeze(img4D(1,:,:,slices(i))),90);
        else
            img(((ceil(i/9)-1)*tamy+1):((ceil(i/9)-1)*tamy+tamy),(mod(i-1,9)*tamx+1):(mod(i-1,9)*tamx+tamx)) = squeeze(img4D(1,:,:,slices(i)));
        end
        
    end
else
    error('Uncorrect dimension order. Please, permute the dimensions until the first one is the singleton');
end