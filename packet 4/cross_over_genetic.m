function [cross_over] = cross_over_genetic(parent,t)
x = parent;
l = size(x);
for i=1:l(1)/2
    x1 = x(i,:);
    x2 = x(l(1)+1-i,:);
    b = randi(l(2)-2)+1;
    y1 = [x1(1:b) x2(b+1:l(2))]; 
    y2 = [x2(1:b) x1(b+1:l(2))]; 
    %mutation
    a = randi(l(2));
    y1(1,a) = mod(randi(t)+y1(1,a),t)+1;
    x(i,:) = y1;
    y2(1,a) = mod(randi(t)+y2(1,a),t)+1;
    x(l(1)+1-i,:) = y2; 
end
cross_over = x;