clear all
clc

% pendulum equation
xn=[0, 0.2, 0.6, 1];
n=10;
h=1/n;
y1(1)=pi/2;
y2(1)=0;

for i=1:n
    y1(i+1) = y1(i)+h*y2(i);
    y2(i+1) = y2(i)-h*sin(y1(i));
end

for j = 2:length(xn)
    ansy1(j-1)=y1(int16(xn(j)/h +1));
    ansy2(j-1)=y2(int16(xn(j)/h +1));
end

%exact solution 
Y1 = @(x) 2*asin((ellipj(ellipticK(1/2)-x,1/2))/sqrt(2));

result=[xn(2:4); ansy1; Y1(xn(2:4)) ;Y1(xn(2:4))-ansy1; ansy2];


%print
h
formatSpec='xn = %.1f \t y_1n = %1.4f \t realY1 = %1.4f \t error_y1=%1.4f \t y_2n = %1.4f \n';
fprintf(formatSpec, result)