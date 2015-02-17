function status = write2file(filename, towrite)
% Writes and appends lines to file. 
if(exist(filename,'file')==2)
    file=fopen(filename,'a');
else
    file=fopen(filename,'w');
end
fprintf(file,'%s\n',towrite);
status=fclose(file);
