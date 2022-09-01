
function [f,df] = surfaceSphere(Nodes,Verts,varargin)
x1 = Verts(2,1); y1 = Verts(2,2);
nNd = size(Nodes,1);
f = zeros(nNd,1);
df = zeros(nNd,3*nNd);
for i = 1:nNd
    x=Nodes(i,2); y=Nodes(i,3);
    [f1,df1] = surfaceFunc(x,y,x1,y1);
%     f(i) = (z-real(f1))^2;
%     df(i,i+nNd*(0:2))=2*(z-real(f1))*([0 0 1]-real(df1));
    f(i) = round(f1,5);
    df(i,i+nNd*(0:2))=df1;
end
if nargin>2
    dNodesdx = varargin{1};
    df = df*reshape(dNodesdx(:,2:4,:),3*size(Nodes,1),size(dNodesdx,3));
end
end

function [f,df] = surfaceFunc(x,y,x1,y1)
% f  = (((x1^2 - 2*x1*x2 + x2^2 + y1^2 - 2*y1*y2 + y2^2)*(x1^2 - 2*x1*x3 + x3^2 + y1^2 - 2*y1*y3 + y3^2)*(x2^2 - 2*x2*x3 + x3^2 + y2^2 - 2*y2*y3 + y3^2))/(4*(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2)^2) - (x - (x1^2*y2 - x1^2*y3 - x2^2*y1 + x2^2*y3 + x3^2*y1 - x3^2*y2 + y1^2*y2 - y1^2*y3 - y1*y2^2 + y1*y3^2 + y2^2*y3 - y2*y3^2)/(2*(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2)))^2 - (y - (- x1^2*x2 + x1^2*x3 + x1*x2^2 - x1*x3^2 + x1*y2^2 - x1*y3^2 - x2^2*x3 + x2*x3^2 - x2*y1^2 + x2*y3^2 + x3*y1^2 - x3*y2^2)/(2*(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2)))^2)^(1/2);
% df = [-(2*x - (x1^2*y2 - x1^2*y3 - x2^2*y1 + x2^2*y3 + x3^2*y1 - x3^2*y2 + y1^2*y2 - y1^2*y3 - y1*y2^2 + y1*y3^2 + y2^2*y3 - y2*y3^2)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2))/(2*(((x1^2 - 2*x1*x2 + x2^2 + y1^2 - 2*y1*y2 + y2^2)*(x1^2 - 2*x1*x3 + x3^2 + y1^2 - 2*y1*y3 + y3^2)*(x2^2 - 2*x2*x3 + x3^2 + y2^2 - 2*y2*y3 + y3^2))/(4*(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2)^2) - (x - (x1^2*y2 - x1^2*y3 - x2^2*y1 + x2^2*y3 + x3^2*y1 - x3^2*y2 + y1^2*y2 - y1^2*y3 - y1*y2^2 + y1*y3^2 + y2^2*y3 - y2*y3^2)/(2*x1*y2 - 2*x2*y1 - 2*x1*y3 + 2*x3*y1 + 2*x2*y3 - 2*x3*y2))^2 - (y - (- x1^2*x2 + x1^2*x3 + x1*x2^2 - x1*x3^2 + x1*y2^2 - x1*y3^2 - x2^2*x3 + x2*x3^2 - x2*y1^2 + x2*y3^2 + x3*y1^2 - x3*y2^2)/(2*x1*y2 - 2*x2*y1 - 2*x1*y3 + 2*x3*y1 + 2*x2*y3 - 2*x3*y2))^2)^(1/2)), -(2*y - (- x1^2*x2 + x1^2*x3 + x1*x2^2 - x1*x3^2 + x1*y2^2 - x1*y3^2 - x2^2*x3 + x2*x3^2 - x2*y1^2 + x2*y3^2 + x3*y1^2 - x3*y2^2)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2))/(2*(((x1^2 - 2*x1*x2 + x2^2 + y1^2 - 2*y1*y2 + y2^2)*(x1^2 - 2*x1*x3 + x3^2 + y1^2 - 2*y1*y3 + y3^2)*(x2^2 - 2*x2*x3 + x3^2 + y2^2 - 2*y2*y3 + y3^2))/(4*(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2)^2) - (x - (x1^2*y2 - x1^2*y3 - x2^2*y1 + x2^2*y3 + x3^2*y1 - x3^2*y2 + y1^2*y2 - y1^2*y3 - y1*y2^2 + y1*y3^2 + y2^2*y3 - y2*y3^2)/(2*x1*y2 - 2*x2*y1 - 2*x1*y3 + 2*x3*y1 + 2*x2*y3 - 2*x3*y2))^2 - (y - (- x1^2*x2 + x1^2*x3 + x1*x2^2 - x1*x3^2 + x1*y2^2 - x1*y3^2 - x2^2*x3 + x2*x3^2 - x2*y1^2 + x2*y3^2 + x3*y1^2 - x3*y2^2)/(2*x1*y2 - 2*x2*y1 - 2*x1*y3 + 2*x3*y1 + 2*x2*y3 - 2*x3*y2))^2)^(1/2)), 0];


zc = -24;
f  = zc + (- x^2 + x1^2 - y^2 + y1^2 + zc^2)^(1/2);
df = [-x/(- x^2 + x1^2 - y^2 + y1^2 + zc^2)^(1/2), -y/(- x^2 + x1^2 - y^2 + y1^2 + zc^2)^(1/2), 0];
end



 
