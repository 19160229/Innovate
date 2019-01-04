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