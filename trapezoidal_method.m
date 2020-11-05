clear all
clc

n=4;
h=1/n;
y1(1)=1; %initial value

y_diff = @(y) -y^2; %diff formulationmtal
Y1 = @(x) 1./(1+x); %exact solution
y_h_2h=[];

for k=1:2;
    h=h*k;
    xn=[0:h:5];
    y1(2)=y1(1)+h*y_diff(y1(1));
        
    for i=1:length(xn)-1
        y1(i+1)=y1(i)+h*y_diff(y1(i));
        y1(i+1)=y1(i)+2*h*y_diff(y1(i+1));
        temp=y1(i+1);
        
        for j=1:100
                y1(i+1)=y1(i)+h*(y_diff(y1(i))+y_diff(temp))/2; %trapezoidal method
            if abs(temp-y1(i+1))<10^-7
                break;
            end
            temp=y1(i+1);
        end
        y1(i+1)=temp;
    end
    
    
    index_x=1;
    for j=1:length(xn)
        if xn(j)==index_x
            y_h_2h=[y_h_2h, y1(j)];
            index_x=index_x+1;
        end
    end
end


y_h_2h=reshape(y_h_2h,[5,2]);

index_x=1;
exact_mat=[];
for i=1:length(xn)
    if xn(i)==index_x
        exact_mat=[exact_mat; xn(i), Y1(xn(i))];
        index_x=index_x+1;
    end
end
result=[exact_mat(:,1), y_h_2h(:,2), exact_mat(:,2)-y_h_2h(:,2), y_h_2h(:,1), exact_mat(:,2)-y_h_2h(:,1),(y_h_2h(:,1)-y_h_2h(:,2))/3];


% error=abs(Y1(xn)-y1);

%print
result=result';
fprintf('x \t\t y_2h \t\t\t Y-y_2h \t\t y_h \t\t\t Y-y_h \t\t\t [y_h-y_2h]/3\n');
fprintf('=======================================================================================\n');
formatSpec='%.1f \t %1.6f \t\t %1.6f \t\t %1.6f \t\t %1.6f \t\t %1.6f\n';
fprintf(formatSpec, result)
