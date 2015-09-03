function [X,w]=gaussquadpts_pickptsol(mu,P,n)
%please keep n(the total number of points in 1D) odd
%mu is column
%this scheme actually solves the set of nonlinear equations by
%picking points on the principle axis of contours that we select


%Dimension of the system
n=5;
N=length(mu);
%% generating the points along the principle axis
x=zeros(n,N);
x(1,:)=[0,0];
[u,s,v]=svd(P);
AA=chol(s)*u;
sqrP=AA';
k=1;
R=[cosd(45),-sind(45);sind(45),cosd(45)];
for i=1:1:(n-1)/2
     if i==1
    for j=1:1:N
    x(k+j,:)=(sqrt(2*i+1)*(1)*sqrP(:,j))';
    x(k+j+N,:)=(-sqrt(2*i+1)*(1)*sqrP(:,j))';
    end
    
    else
        
    %for j=1:1:N
    x(k+1,:)=sqrt(2*i+1)*(sqrP(:,1)+sqrP(:,2))';
    x(k+2,:)=sqrt(2*i+1)*(sqrP(:,1)-sqrP(:,2))';
    x(k+3,:)=sqrt(2*i+1)*(-sqrP(:,1)-sqrP(:,2))';
    x(k+4,:)=sqrt(2*i+1)*(-sqrP(:,1)+sqrP(:,2))';
   % end
    end
    
    k=k+2*N;
end
% size(x)
plot(x(:,1),x(:,2),'bo')
hold on
plot([sqrP(1,1),-sqrP(1,1)],[sqrP(1,2),-sqrP(1,2)])
plot([sqrP(1,2),-sqrP(1,2)],[sqrP(2,2),-sqrP(2,2)])
axis([-15,15,-15,15])
%% generating moments of the given PDF
% only generating even moments as odd moments are already 
% satisfied by symmetric points
%1+(n-1)/2  weights are there.
%(n-1)/2 even moments are required as sum(w)=1 is another constraint

m=[1;P(1,1)];%P(1,2);P(2,2)];
mm=permute_moments(P,4);
m=[m;mm(1)];
    
%% now generating the A matrix of Aw=m equation.
% A=ones(1,3);
A=[1,4,4];
%second order moment
for i=0:1:0
    a=(x(:,1).^(2-i)).*(x(:,2).^(i));
    a=a';
    b=[0];
    for k=2:2*N:length(x)
        b=[b,sum(a(k:1:k+2*N-1))];
    end
 
    A=vertcat(A,b);
end
%fourth order moment
for i=0:1:0
    a=(x(:,1).^(4-i)).*(x(:,2).^(i));
   
    a=a';
    b=[0];
    for k=2:2*N:length(x)
        b=[b,sum(a(k:1:k+2*N-1))];
    end
    A=vertcat(A,b);
  
end
%  eig(A);
%% calculating the weights of the system
w=A\m;
% sum(w)
% m
% w=0;
% [A(3,:)./A(2,:)]
% [A(3,:)./A(4,:)]
% 
% [A(5,:)./A(6,:)]
% [A(6,:)./A(7,:)]
% [A(7,:)./A(8,:)]
% [A(8,:)./A(9,:)]
% [A(7,:)./A(9,:)]
% [A(7,:)./A(5,:)]
% [A(8,:)./A(5,:)]
% [A(9,:)./A(5,:)]
% [A(6,:)./A(8,:)]
% [A(6,:)./A(9,:)]
% m
% A
% [m(3)/m(2),m(3)/m(4),m(7)/m(8)]
%% now shfit the PDF to the required mean.
X=x;
w=[w(1),w(2)*ones(1,4),w(3)*ones(1,4)];
w=w';
% A
m
[sum(w.*(X(:,1).^(2)).*(X(:,2).^(0)));...
sum(w.*(X(:,1).^(1)).*(X(:,2).^(1)));...
sum(w.*(X(:,1).^(0)).*(X(:,2).^(2)));...
sum(w.*(X(:,1).^(4)).*(X(:,2).^(0)));...
sum(w.*(X(:,1).^(3)).*(X(:,2).^(1)));...
sum(w.*(X(:,1).^(2)).*(X(:,2).^(2)));...
sum(w.*(X(:,1).^(1)).*(X(:,2).^(3)));...
sum(w.*(X(:,1).^(0)).*(X(:,2).^(4)))]
%X=x+[ones(mu(1)*length(x),1),mu(2)*ones(length(x),1)];