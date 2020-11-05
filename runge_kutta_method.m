clear all
clc


y_1 = @(x,y) 1/(1+x^2) - 2*y^2;

Y = @(x) x/(1+x^2);

h0=1/2;
err=[];
for i=1:5
    y0 = 0;
    h(i)=h0^(i+2);
    x=[0:h(i):2];
    fprintf('\n\n')
    fprintf('i = %1.0f \t h = %1.5f\n', i+2, h(i));
    fprintf('=============================================\n');
    temp=0;
    for j=1:length(x)-1;
        v1=y_1(x(j),y0);
        v2=y_1(x(j)+h(i)/2,y0+h(i)*v1/2);
        v3=y_1(x(j)+h(i)/2,y0+h(i)*v2/2);
        v4=y_1(x(j)+h(i),y0+h(i)*v3);
        y0=y0+h(i)*(v1+2*v2+2*v3+v4)/6;
        
        error=Y(x(j+1))-y0;
        if mod(x(j+1),2)==0 && x(j)>0
            temp=temp+1;
            fprintf('x = %1.0f \t y = %1.8f \t error = %1.1e \n',x(j),y0,error);
            err(temp,i)=error;
        end
        
        
        
    end
    
end
for k=1:4;
    ratio=err(:,k)./err(:,k+1);
    h_2h=(err(:,k)-err(:,k+1))/15;
    fprintf('\n\n')
    fprintf('h = %1.5f \t 2h = %1.5f\n', (1/2)^(k+3),(1/2)^(k+2))
    fprintf('------------------------------\n')
    fprintf('ratio = %2.2f \n',ratio)
    fprintf('(y_h - y_2h)/15 = %1.1e \n',h_2h)
end