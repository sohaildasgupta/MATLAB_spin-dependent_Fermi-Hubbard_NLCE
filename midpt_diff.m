function array = midpt_diff(x,y)
a = [];
for i=2:numel(y)-1
    a(end+1) = (x(i)-x(i+1))/((x(i-1)-x(i))*(x(i-1)-x(i+1)))*y(i-1)...
    +(2*x(i)-x(i-1)-x(i+1))/((x(i)-x(i-1))*(x(i)-x(i+1)))*y(i)...
    +(x(i)-x(i-1))/((x(i+1)-x(i-1))*(x(i+1)-x(i)))*y(i+1);
end
array = a;
end