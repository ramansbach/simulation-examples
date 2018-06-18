function [fig] = errorMesh(x,y,z,zerrs,color,fig,dx,emptyinds)
%function that essentially does surf except it also shows errors in the z
%direction and it takes 1D vectors
xu = unique(x);
yu = unique(y);

view(3);
grid on;
for i = 1:length(xu)
    xi = find(x==xu(i));
    plot3(x(xi),y(xi),z(xi),'-','linewidth',2,'color',color)
end

for i = 1:length(yu)
   yi = find(y==yu(i));
   plot3(x(yi),y(yi),z(yi),'-','linewidth',2,'color',color);
end

for i = 1:length(x)
   plot3([x(i) x(i)],[y(i) y(i)],[z(i)-zerrs(i) z(i)+zerrs(i)],'-','color',color,'linewidth',1.5) 
end

for i = 1:length(x)
   plot3([x(i)-dx x(i) + dx],[y(i) y(i)],[z(i)-zerrs(i) z(i) - zerrs(i)],'-','color',color,'linewidth',1.5);
   plot3([x(i)-dx x(i) + dx],[y(i) y(i)],[z(i)+zerrs(i) z(i) + zerrs(i)],'-','color',color,'linewidth',1.5);
end

if exist('emptyinds','var')
   scatter3(x(emptyinds{1}),y(emptyinds{1}),z(emptyinds{1}),100,z(emptyinds{1}),'linewidth',3);
   scatter3(x(emptyinds{2}),y(emptyinds{2}),z(emptyinds{2}),100,z(emptyinds{2}),'filled');
else
    scatter3(x,y,z,100,z,'filled')
end