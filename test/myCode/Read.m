clc
clear

%����Ҷ��ļ���
directoryM = uigetdir('','ѡ��ָ���ļ��У�');
dirsM=dir(directoryM);%dirs�ṹ������,���������ļ������������ļ�������Ϣ��
dircellM=struct2cell(dirsM)'; %����ת����ת��ΪԪ������
filenamesM=dircellM(:,1) ;%�ļ����ʹ���ڵ�һ��
%Ȼ����ݺ�׺��ɸѡ��ָ�������ļ�������

%�����Ե�ļ���
directoryN = uigetdir('','ѡ��ָ���ļ��У�');
dirsN=dir(directoryN);%dirs�ṹ������,���������ļ������������ļ�������Ϣ��
dircellN=struct2cell(dirsN)'; %����ת����ת��ΪԪ������
filenamesN=dircellN(:,1) ;%�ļ����ʹ���ڵ�һ��

[n m] = size(filenamesM);%��ô�С

cnt=0;
for i = 1:n
    if ~isempty( strfind(filenamesM{i}, '.csv') ) && ~isempty( strfind(filenamesN{i}, '.csv') )%ɸѡ��csv�ļ�
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
%d��ʾ��Ҫ������ά����k��ʾ���ڸ�����Ҫ��d<k����kҪС��һ����Ƶ����֡����
%���d<=k-2��������warning��
[T,NI] = LITSA(res',6,8);
T=T';
save([directoryM '.mat'],'T');