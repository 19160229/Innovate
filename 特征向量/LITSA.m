function [T,NI] = LITSA(data,d,K,NI)%增量特征提取
%Generate sampled data[2;2;3]

[m,N] = size(data);  % m is the dimensionality of the input sample points. 
% Step 0:  Neighborhood Index
if nargin<4
    if length(K)==1
        K = repmat(K,[1,N]);
    end;
    NI = cell(1,N); 
    if m>N
        a = sum(data.*data); 
        dist2 = sqrt(repmat(a',[1 N]) + repmat(a,[N 1]) - 2*(data'*data));
        for i=1:N
            % Determine ki nearest neighbors of x_j
            [dist_sort,J] = sort(dist2(:,i));  
            Ii = J(1:K(i)); 
            NI{i} = Ii;
        end;
    else
        for i=1:N
            % Determine ki nearest neighbors of x_j
            x = data(:,i); ki = K(i);
            dist2 = sum((data-repmat(x,[1 N])).^2,1);    
            [dist_sort,J] = sort(dist2);  
            Ii = J(1:ki);  
            NI{i} = Ii;
        end;
    end;
else
    K = zeros(1,N);
    for i=1:N
        K(i) = length(NI{i});
    end;
end;
% Step 1:  local information:ipca
BI={}; 
Thera = {}; 

%  load('mnist_all');   
%  %a(1:10,1:784)=train0(1:10,:);
%  a(1:10,1:784)=test2(1:10,:);
%  a=a(1,:);
%  a=a';
%  size(a)
 
% load('PIE_32x32.mat'); 
% load('fea'); 
% a=fea(1,:)';
% size(a)

% load('ar20_50x40.mat'); 
% %load('fea'); 
% a=fea(1,:)';
% size(a)

  load('Yale_32x32.mat');
  a=fea(1,:)';
  size(a)
    
%   load('Yale_64x64.mat');
%   a=fea(1,:)';
%   size(a)

%   load('YaleB_32x32.mat');
%   a=fea(1,:)';
%   size(a)

%  load('usps_all');
%  a(1:256,1:1100)=data(:,:,1);
%  a=a(:,1);
%  a=double(a);
%  size(a)
%  a=a';

%  load('face_data'); 
%  a=images(:,1);
%  a=images(1,:)';
 
%   load('UMIST.mat');
%   a=fea(1,:)';
%   size(a)
 
%    load('frey_rawface'); 
%    [D,N]=size(ff);
%    N=600;
%    fff=zeros(D,N);
%     for i=1:D
%         for j=601:N+600
%             fff(i,j-600)=ff(i,j-600);
%             if fff(i,j-600)>127
%             fff(i,j-600)=fff(i,j-600)-256;
%             end
%         end
%     end
%     fff=fff-repmat(sum(fff,2)/N,1,N);
%     a=fff(:,1);
%     size(fff)   
    
    
for i=1:N
    % Compute the d largest right singular eigenvectors of the centered matrix
    Ii = NI{i};
    ki = K(i);
    M=repmat(mean(data(:,Ii),2),[1,ki]);
    t=(N-1)/N;
    %M=t*M+(1-t)*repmat([2;2;3],[1,ki]);
    %newdata=load('newface.m');%新加入点的坐标
    newdata=a;
    %a=a';
    
    Xi=data(:,Ii)-repmat(mean(data(:,Ii),2),[1,ki]);
    C=Xi'*Xi;
    
    %C_n=repmat([2;2;3],[1,ki])-M;
    C_n=repmat(newdata,[1,ki])-M;
    
    C=t*C+(1-t)*(C_n'*C_n);
    opts.disp=0;
    [Vi,Si]=eigs(C,d+1,'lm',opts); 
    %2~d+1??
    Vi=Vi(:,1:d+1);
    
    Li=diag(Si);
    Ai=[Vi*sqrt(t*Li),sqrt(1-t).*(C_n')];
%end;

B=Ai*Ai';
B=(B+B')/2;
[Ui,Di]=schur(B);
[s,Ji]=sort(diag(Di));
Ui=Ui(:,Ji(1:d));
[m,n]=size(Ui);
    Gi=[repmat(1/sqrt(ki),[m,1]) Ui];
    BI{i}=eye(ki)-Gi*Gi';
end;
D=speye(N);
for i=1:N
    Ii=NI{i};
    D(Ii,Ii)=D(Ii,Ii)+BI{i};
    D(i,i)=D(i,i)-1;
end;

D=(D+D')/2;
options.disp=0;
options.isreal=1;
options.issym=1;
[Y,E]=eigs(D,d+2,0,options);
lambda=diag(E);
[lambda_s,J]=sort(abs(lambda));
Y=Y(:,J);lambda=lambda(J);

%d=2;K=8;m=3;t=(N-1)/N;
T=Y(:,2:d+1)';
save Y1.m T -ASCII;