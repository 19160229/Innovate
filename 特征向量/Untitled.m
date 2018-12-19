clear;
m=[1	-1	1	-1	-1	1	-1	-1];
len=1208;
for i=1:len
    Y(i,:)=m;
end
save('Y.mat','Y');