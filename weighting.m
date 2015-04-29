function y = weighting(p,d,k)
sumatoria = 1;
for jarl=1:length(d),
    if ndims(p)==1
        sumatoria = sumatoria + power(-1,1-d(jarl))*exp(-squeeze(abs(p(jarl)))/k);
    elseif ndims(p)==2
        sumatoria = sumatoria + power(-1,1-d(jarl))*exp(-squeeze(abs(p(jarl,:)))/k);
    elseif ndims(p)==3
        sumatoria = sumatoria + power(-1,1-d(jarl))*exp(-squeeze(abs(p(jarl,:,:)))/k);
    else
        error('P matrices with dimensions > 3 not supported');
    end
end
y = sumatoria/2; 
end
