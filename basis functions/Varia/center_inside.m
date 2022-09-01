function [in] = center_inside(Nodes,Trigs)
C = @(t,c) Nodes(Trigs(:,t+1),c+1);
x1=C(1,1); y1=C(1,2); x2=C(2,1); y2=C(2,2); x3=C(3,1); y3=C(3,2);

xc=(x1.^2.*y2-x1.^2.*y3-x2.^2.*y1+x2.^2.*y3+x3.^2.*y1-x3.^2.*y2+y1.^2.*y2-y1.^2.*y3-y1.*y2.^2+y1.*y3.^2+y2.^2.*y3-y2.*y3.^2)./(2.*(x1.*y2-x2.*y1-x1.*y3+x3.*y1+x2.*y3-x3.*y2));
yc=(-x1.^2.*x2+x1.^2.*x3+x1.*x2.^2-x1.*x3.^2+x1.*y2.^2-x1.*y3.^2-x2.^2.*x3+x2.*x3.^2-x2.*y1.^2+x2.*y3.^2+x3.*y1.^2-x3.*y2.^2)./(2.*(x1.*y2-x2.*y1-x1.*y3+x3.*y1+x2.*y3-x3.*y2));

k = convhull(Nodes(:,2),Nodes(:,3));
in = inpolygon(xc,yc,Nodes(k,1),Nodes(k,2));
end