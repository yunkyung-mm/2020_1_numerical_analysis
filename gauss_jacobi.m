clear all
clc

A=[10 3 1; 2 -10 3; 1 3 10];
b=[14;-5;14];
x=[0;0;0];
new_A=A;

iter=0;
error(1)=max(abs([1;1;1]-x));
fprintf('iter = %d, x(1) = %0.6f, x(2) = %0.6f, x(3) = %0.6f, error = %0.2f\n',iter, x', error(1))

diag_A=[10;-10;10];
for iter=1:6
    for i=1:length(x);       
        new_A(i,i)=0;
        x(i)=(b(i)-new_A(i,:)*x)./diag_A(i);
    end
    error(iter+1)=max(abs([1;1;1]-x));
    ratio(iter)=error(iter+1)/error(iter);
    fprintf('iter = %d, x(1) = %0.6f, x(2) = %0.6f, x(3) = %0.6f, error = %1.2e, ratio = %0.2f \n',iter, x', error(iter+1),ratio(iter))
    
end