% clear all
% clc

n=16;
h=1/n;
xn=[0:h:2.25];
lambda=1;
y1(1)=1; %initial value
y1(2)=y1(1)+lambda*h*y1(1); %Euler method

for i=2:length(xn)-1
    y1(i+1) = y1(i-1)+2*h*lambda*y1(i); %mid point method 
end

%exact solution 
Y1 = @(x) exp(lambda*x);

result=[xn; y1; Y1(xn) ;Y1(xn)-y1];
error=abs(Y1(xn)-y1);

%print
h
formatSpec='xn = %.2f \t\t y_n = %1.4f \t\t realY = %1.4f \t\t error=%1.4f \n';
fprintf(formatSpec, result)



plot(xn, log2(error))
hold on 
% plot(xn, Y1(xn))
legend('n=4','n=8', 'n=16')
title('log2(error) graph')
% hold off