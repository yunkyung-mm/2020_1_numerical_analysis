clear all
clc

h=[0.125, 0.0625];
y_diff=@(x,y) 1/(1+x^2) - 2*y^2;
y_exact=@(x) x/(1+x^2);

y(1) = 0;

for k=1:2
    x=[0:h(k):10];
    n=length(x);
    y(2)=y_exact(x(2));
    y(3)=y_exact(x(3));
    fprintf('h = %1.4f \n',h(k))
    fprintf('------------------\n')
    temp=1;
    for i = 3 : n % main phase
        y(i+1) = y(i)+(23*y_diff(x(i),y(i))-16*y_diff(x(i-1),y(i-1))+5*y_diff(x(i-2),y(i-2)))*h(k)/12;
        x(i+1) = x(i)+h(k);
        y(i+1) = y(i)+(9*y_diff(x(i+1),y(i+1))+19*y_diff(x(i),y(i))-5*y_diff(x(i-1),y(i-1))+y_diff(x(i-2),y(i-2)))*h(k)/24;
        error=y_exact(x(i+1))-y(i+1);
        
        if mod(x(i+1),2)==0 && x(i)>0;
            errorx(k,temp)=error;
            fprintf('%1.0f\t%1.2e\n', x(i+1), error)
            temp=temp+1;
        end
    end
    fprintf('\n\n')
end

fprintf('ratio\n')
fprintf('------------------\n')
for i=1:length(errorx)
    fprintf('%2.1f\n',errorx(1,i)/errorx(2,i))
end