function outputMatrix(x, d1, d2)
fs = strcat(' %',num2str(d1),'.', num2str(d2),'f');
fs1 = strcat(fs,' &  ');
fs2 = strcat(fs,' \\\\ \n');

for ir = 1:size(x,1)
    fprintf(1,fs1,x(ir,1:end-1));
    fprintf(1,fs2,x(ir,end));
end