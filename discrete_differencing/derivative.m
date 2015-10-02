function [D]=derivative(f)

%%finds the 5 point stencil central differences first derivative for
%%equally spaced samples

for i=3:length(f)-3;
    D(i)=(-f(i+2)+8*f(i+1)-8*f(i-1)+f(i-2))/12;
end

for i=2;D(i)=(f(i)+f(i+1))/2-(f(i-1)+f(i))/2; end
 
for i=length(f)-2; D(i)=(f(i)+f(i+1))/2-(f(i-1)+f(i))/2; end

for i=1; D(i)=f(i+1)-f(i);end
for i=length(f)-1; D(i)=f(i+1)-f(i);end

