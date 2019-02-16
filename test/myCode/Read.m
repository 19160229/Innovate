clc
clear

%读入灰度文件夹
directoryM = uigetdir('','选择指定文件夹：');
dirsM=dir(directoryM);%dirs结构体类型,不仅包括文件名，还包含文件其他信息。
dircellM=struct2cell(dirsM)'; %类型转化，转化为元组类型
filenamesM=dircellM(:,1) ;%文件类型存放在第一列
%然后根据后缀名筛选出指定类型文件并读入

%读入边缘文件夹
directoryN = uigetdir('','选择指定文件夹：');
dirsN=dir(directoryN);%dirs结构体类型,不仅包括文件名，还包含文件其他信息。
dircellN=struct2cell(dirsN)'; %类型转化，转化为元组类型
filenamesN=dircellN(:,1) ;%文件类型存放在第一列

[n m] = size(filenamesM);%获得大小

cnt=0;
for i = 1:n
    if ~isempty( strfind(filenamesM{i}, '.csv') ) && ~isempty( strfind(filenamesN{i}, '.csv') )%筛选出csv文件
        cnt=cnt+1;
        filenameM = filenamesM{i};
        filepathM = fullfile(directoryM,filenameM);
        filenameN = filenamesN{i};
        filepathN = fullfile(directoryN,filenameN);
        tmpM=csvread(filepathM);
        tmpM(:,1281)=[];
        tmpN=csvread(filepathN);
        tmpN(:,1281)=[];
        tmpN((isnan(tmpN)==1)) = 0;
%         M(:,:,cnt)=tmp;
        res(cnt,:)=[reshape(tmpM,[],1)',reshape(tmpN,[],1)'];
 %       disp([filepath(1:end-4) '.mat']);
    end
end
clear tmpM;
clear tmpN;

%LITSA(res',d,k)    
%d表示想要降到的维数，k表示近邻个数，要求d<k，而k要小于一个视频的总帧数，
%最好d<=k-2，否则有warning。
[T,NI] = LITSA(res',6,8);
T=T';
save([directoryM '.mat'],'T');