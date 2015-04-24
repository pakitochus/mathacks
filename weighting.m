function y = weighting(p1,p2,k)
y = ((exp(-k*(p1-0.5))-exp(k*((1-p2)-0.5)))/exp(k/2)+1)/2;
end