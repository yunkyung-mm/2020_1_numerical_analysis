function detrap(f, x_0, y_0, x_end, eps, h_init, h_min, h_max, ier)
% 1 initial value setting 
Y=@(x) x/(1+x^2);
y_h(1)=y_0;
y_2h(1)=y_0;
h(1)=h_init;
x(1)=x_0;

% 2 initialize
loop=1;
ier=0;

% 3 choose an initial value of h 
i=1;
fprintf('x_n \t h \t\t y_n\t\tY(x_n)-y_n \ttrunc \n')

while 1
    % 4 calculate y_h(x_0+h) y_h/2(x_0+(h/2)) y_h/2(x_0+h)
    % using 1 euler method and 2 trapezoidal method 
    y_hp=y_h(i)+h(i)*f(x(i),y_h(i)); %euler 
    y_hp=y_h(i)+h(i)*(f(x(i),y_h(i))+f(x(i)+h(i),y_hp))/2; %trapezoidal
    y_h(i+1)=y_h(i)+h(i)*(f(x(i),y_h(i))+f(x(i)+h(i),y_hp))/2; %trapezoidal
    
    temp2=y_2h(2*i-1)+h(i)*f(x(i),y_2h(2*i-1))/2; %euler
    temp2=y_2h(2*i-1)+h(i)*(f(x(i),y_2h(2*i-1))+f(x(i)+h(i)/2,temp2))/4; %trapezoidal
    y_2h(2*i)=y_2h(2*i-1)+h(i)*(f(x(i),y_2h(2*i-1))+f(x(i)+h(i)/2,temp2))/4; %trapezoidal
    
    temp3=y_2h(2*i)+h(i)*f(x(i)+h(i)/2,y_2h(2*i))/2; %euler
    temp3=y_2h(2*i)+h(i)*(f(x(i)+h(i)/2,y_2h(2*i))+f(x(i)+h(i),temp3))/4; %trapezoidal
    y_2h(2*i+1)=y_2h(2*i)+h(i)*(f(x(i)+h(i)/2,y_2h(2*i))+f(x(i)+h(i),temp3))/4; %trapezoidal
    
    % 5 calculate truncation error
    trunc(i)=4/3*(y_2h(2*i+1)-y_h(i+1));
    
    % 6 
    if (eps*h(i)/4<=abs(trunc(i))) && (abs(trunc(i)) <=eps*h(i))
        loop=2;
        h(i+1)=h(i);
        x(i+1)=x(i)+h(i);
        fprintf('%0.4f\t%0.4f\t%0.6f\t%1.2e\t%1.2e\n',x(i+1),h(i+1),y_h(i+1),Y(x(i+1))-y_h(i+1),trunc(i))
        break
    else
        % 7 calculate D_3y
        f0=f(x(i),y_h(i));
        f1=f(x(i)+h(i)/2,y_2h(2*i));
        f2=f(x(i)+h(i),y_2h(2*i+1));
        d3y=(f2-2*f1+f0)/(h(i)/2)^2;
        
        if d3y==0
            loop=2;
            h(i)=h_max;
        else
            h(i)=sqrt(6*eps/abs(d3y));
        end
    end
    
end


i=2;
% 10 predictor-corrector step 
while 1
    while 1
        % 11 
        temp=y_h(i-1)+2*h(i)*f(x(i),y_h(i)); %midpoint
        y_h(i+1)=y_h(i)+(h(i)/2)*(f(x(i),y_h(i)) + f(x(i)+h(i),temp));%trapezoidal
        
        % 12 
        trunc(i)=-(y_h(i+1)-temp)/6;
      
        % 13
        if abs(trunc(i))<eps*h(i)/4 || eps*h(i)<abs(trunc(i)) 
            break
        end
        
        h(i+1)=h(i);
        x(i+1)=x(i)+h(i);
        
        % 14 print 
        fprintf('%0.4f\t%0.4f\t%0.6f\t%1.2e\t%1.2e\n',x(i+1),h(i),y_h(i+1),Y(x(i+1))-y_h(i+1),trunc(i))
        i=i+1;
    end
    
    % 16 change the stepsize
    
    % 17 
    h_new=sqrt(eps*h(i)^3/(2*abs(trunc(i))));
    
    % 18 
    h(i)=min(h_new,2*h(i));
    
    % 19 
    if h(i)<h_min
        ier=2;
        break
    elseif h(i)>h_max
        ier=1;
        h(i)=h_max;
    end
    
    % 20 
    y0=y_h(i)+h(i)*f(x(i),y_h(i)); %euler
    y0=y_h(i)+h(i)*(f(x(i),y_h(i))+f(x(i)+h(i),y0))/2; %trapezoidal
    y_h(i+1)=y_h(i)+h(i)*(f(x(i),y_h(i))+f(x(i)+h(i),y0))/2;%trapezoidal
    
    h(i+1)=h(i);
    x(i+1)=x(i)+h(i);
    
    % 21 print
    fprintf('%0.4f\t%0.4f\t%0.6f\t%1.2e\t%1.2e\n',x(i+1),h(i+1),y_h(i+1),Y(x(i+1))-y_h(i+1),trunc(i))   
    i=i+1;
    
    % 22 
    if x(i)>x_end
        break
    end
    
end


end





