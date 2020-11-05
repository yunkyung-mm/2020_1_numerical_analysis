clear all
clc

u=@(x,t) exp(-0.1*t)*sin(pi*x); %exact
d_0=0;
d_1=0;
f=@(x) sin(pi*x);
G=@(x,t) (pi^2-0.1)*exp(-0.1*t)*sin(pi*x);

m=[4,8,16];


for i=1:3;
    temp=1;
    fprintf('\n\n')
    fprintf('m = %0.1f\n',m(i));
    fprintf('=================\n')
    del=1/m(i);%space size
    h=0.1; %time size
    x=[0:del:1];
    t=[0:h:5];
    u_0=f(x(2:length(x)-1));
    
    lambda=-2*eye(m(i)-1);
    for j=1:m(i)-2;
        lambda(j,j+1)=1;
        lambda(j+1,j)=1;
    end
    lambda=(1/del^2)*lambda;
    V=u_0';
    for k=2:length(t);
        g=G(x(2:length(x)-1),t(k)); 
        V=inv(eye(m(i)-1)-h*lambda)*(V+h*g');
        exact=u(x(2:length(x)-1),t(k));
        
        error=max(abs(exact'-V));
        if t(k)>0 && mod(t(k),1)==0;
            E(temp,i)=error;
            fprintf('t=%0.f \t %0.2e\n',t(k),error)
            temp=temp+1;
        end
    end
    
end

ratio1=E(:,1)./E(:,2);
ratio2=E(:,2)./E(:,3);




