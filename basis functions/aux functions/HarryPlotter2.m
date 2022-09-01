function HarryPlotter2(Nodes,Elements,Sections,varargin)
Bars = Elements(:,[1,5,6]);
xSec = 50*Sections(Elements(:,3),2);
Node = Nodes(setdiff(Nodes(:,1),unique(Elements(:,7))),:);
Tri = delaunay(Node(:,2),Node(:,3));
Tri = [Node(Tri(:,1),1),Node(Tri(:,2),1),Node(Tri(:,3),1)];
figure
trisurf(Tri, Nodes(:,2), Nodes(:,3), Nodes(:,4), 'visible', 'off','facecolor',[0.5,0.5,0.5]);
axis equal; colormap winter; axis off; 
% view(45,35
view(2)
if nargin > 4
    TitleText = varargin{1};
    title(TitleText)
end
ax = gca;
ax.TitleHorizontalAlignment = 'left';
hold on
for i = 1:size(Bars,1)
    x1 = Nodes(Bars(i,2),2); x2 = Nodes(Bars(i,3),2);
    y1 = Nodes(Bars(i,2),3); y2 = Nodes(Bars(i,3),3);
    z1 = Nodes(Bars(i,2),4); z2 = Nodes(Bars(i,3),4);
   line([x1,x2],[y1,y2],[z1,z2],'LineWidth',xSec(i),'color',[0,0,0]); 
end
hold off
drawnow
end