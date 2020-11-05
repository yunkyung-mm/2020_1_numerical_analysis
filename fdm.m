clear all
clc
%u_xx+u_yy=0
%u(0,y)=cos(pi*y)
%u(x,0)=exp(pi*x)
%u(1,y)=exp(pi)*cos(pi*y)
%u(x,1)=-exp(pi*x)
%u(x,y)=exp(pi*x)*cos(pi*y)

u_0y = @(y) cos(pi*y);
u_x0 = @(x) exp(pi*x);
u_1y = @(y) exp(pi)*cos(pi*y);
u_x1 = @(x) -exp(pi*x);
u = @(x,y) exp(pi*x)*cos(pi*y);

u_h = @(x,y) (1-x)*u_0y(y)+x*u_1y(y)+(1-y)*u_x0(x)+y*u_x1(x)+...
    -((1-y)*(1-x)+(1-y)*x*exp(pi)+y*(1-x)*-1 + x*y*-exp(pi));


N=32;
h=1/N;
u_0=zeros(N+1,N+1);
u_0(1:N+1,1)=u_x0(0:h:1);
u_0(1,1:N+1)=u_0y(0:h:1);
u_0(end,1:N+1)=u_1y(0:h:1);
u_0(1:N+1,end)=u_x1(0:h:1);
u_exact=u_0;
for j =2:N;
    for k = 2:N;
        u_0(j,k)=u_h((j-1)*h,(k-1)*h);
    end
end



U(:,:,1)=u_0;
fprintf('SOR method\n')
epsilon=[0.01, 0.001, 0.0001];
for i=1:length(epsilon)
    error=1;
    m=2;
    z=1-2*(sin(pi/(2*N)))^2;
    om=2/(1+sqrt(1-z^2));
    while error>epsilon(i)
        U(1:N+1,1,m)=u_x0(0:h:1);
        U(1,1:N+1,m)=u_0y(0:h:1);
        U(N+1,1:N+1,m)=u_1y(0:h:1);
        U(1:N+1,N+1,m)=u_x1(0:h:1);
        
        for j=2:N
            for k=2:N
                v=1/4*( U(j+1,k,m-1)+ U(j,k+1,m-1)+U(j-1,k,m)+U(j,k-1,m));
                U(j,k,m)=om*v+(1-om)*U(j,k,m-1);
                u_exact(j,k)=u((j-1)*h,(k-1)*h);  
            end
        end
        if m>3
            c=max(max(abs(U(:,:,m)-U(:,:,m-1))))/max(max(abs(U(:,:,m-1)-U(:,:,m-2))));
            error=(c/(1-c))*max(max(abs(U(:,:,m)-U(:,:,m-1))));            
        end
        m=m+1;
    end
    fprintf('N = %d, eps = %0.4f, m = %d \n', N, epsilon(i), m-2)
end
