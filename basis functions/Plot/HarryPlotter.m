function [] = HarryPlotter(Nodes, Trigs, Bars, xSec)
figure
trisurf(Trigs, Nodes(:,2), Nodes(:,3), Nodes(:,4), 'visible', 'on');
axis equal; colormap winter; axis off;
hold on
for i = 1:size(Bars,1)
    x1 = Nodes(Bars(i,1),2); x2 = Nodes(Bars(i,2),2);
    y1 = Nodes(Bars(i,1),3); y2 = Nodes(Bars(i,2),3);
    z1 = Nodes(Bars(i,1),4); z2 = Nodes(Bars(i,2),4);
   line([x1,x2],[y1,y2],[z1,z2],'Color',[1 1 1]*(1-xSec(i)));
end
end

