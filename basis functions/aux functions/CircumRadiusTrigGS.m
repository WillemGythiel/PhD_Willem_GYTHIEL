function [Re] = CircumRadiusTrigGS(Nodes,Trigs,varargin)
nTrig = size(Trigs,1);
if nargin < 3
    nC = 1;
    dNd = 1;
else
    nC = 3*size(Nodes,1);
    dNodesdx = varargin{1};
    dNd = reshape(dNodesdx(:,2:4,:),3*size(Nodes,1),size(dNodesdx,3));
end
Av = zeros(nTrig,nC);
for i = 1:nTrig
    NdT  = Trigs(i,2:4);
    Cor  = Nodes(NdT,2:4);
    x1 = Cor(1); y1 = Cor(4); z1 = Cor(7);
    x2 = Cor(2); y2 = Cor(5); z2 = Cor(8);
    x3 = Cor(3); y3 = Cor(6); z3 = Cor(9);
    if nargin < 3
        j = 1;
        R = (-((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)^(-1.0)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)^(-1.0)*(-((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)+((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)+((x2-x3)^2+(y2-y3)^2+(z2-z3)^2))^2+4)^(-0.5)*((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)^0.5;
    else
        j = Addkron(NdT,size(Nodes,1)*(0:2));
        dfdx1 = ((4*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)-(-((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)+((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)+((x2-x3)^2+(y2-y3)^2+(z2-z3)^2))^2)*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)^(-1.0)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)^(-1.0))^(-2.0)*(((4*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)-(-((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)+((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)+((x2-x3)^2+(y2-y3)^2+(z2-z3)^2))^2)*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)^(-1.0)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)^(-1.0))^0.5*((-x1+x3)*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)*(-((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)+((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)+((x2-x3)^2+(y2-y3)^2+(z2-z3)^2))+2.0*(x2-x3)*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)^2.0)*((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)*(-((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)+((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)+((x2-x3)^2+(y2-y3)^2+(z2-z3)^2))+((4*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)-(-((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)+((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)+((x2-x3)^2+(y2-y3)^2+(z2-z3)^2))^2)*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)^(-1.0)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)^(-1.0))^1.5*(x1-x2)*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)^3.0*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)^2.0)*((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)^(-0.5)*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)^(-3.0)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)^(-2.0);
        dfdy1 = ((4*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)-(-((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)+((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)+((x2-x3)^2+(y2-y3)^2+(z2-z3)^2))^2)*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)^(-1.0)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)^(-1.0))^(-2.0)*(((4*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)-(-((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)+((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)+((x2-x3)^2+(y2-y3)^2+(z2-z3)^2))^2)*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)^(-1.0)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)^(-1.0))^0.5*((-y1+y3)*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)*(-((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)+((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)+((x2-x3)^2+(y2-y3)^2+(z2-z3)^2))+2.0*(y2-y3)*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)^2.0)*((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)*(-((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)+((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)+((x2-x3)^2+(y2-y3)^2+(z2-z3)^2))+((4*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)-(-((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)+((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)+((x2-x3)^2+(y2-y3)^2+(z2-z3)^2))^2)*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)^(-1.0)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)^(-1.0))^1.5*(y1-y2)*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)^3.0*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)^2.0)*((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)^(-0.5)*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)^(-3.0)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)^(-2.0);
        dfdz1 = ((4*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)-(-((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)+((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)+((x2-x3)^2+(y2-y3)^2+(z2-z3)^2))^2)*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)^(-1.0)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)^(-1.0))^(-2.0)*(((4*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)-(-((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)+((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)+((x2-x3)^2+(y2-y3)^2+(z2-z3)^2))^2)*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)^(-1.0)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)^(-1.0))^0.5*((-z1+z3)*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)*(-((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)+((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)+((x2-x3)^2+(y2-y3)^2+(z2-z3)^2))+2.0*(z2-z3)*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)^2.0)*((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)*(-((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)+((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)+((x2-x3)^2+(y2-y3)^2+(z2-z3)^2))+((4*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)-(-((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)+((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)+((x2-x3)^2+(y2-y3)^2+(z2-z3)^2))^2)*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)^(-1.0)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)^(-1.0))^1.5*(z1-z2)*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)^3.0*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)^2.0)*((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)^(-0.5)*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)^(-3.0)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)^(-2.0);
        dfdx2 = ((4*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)-(-((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)+((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)+((x2-x3)^2+(y2-y3)^2+(z2-z3)^2))^2)*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)^(-1.0)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)^(-1.0))^(-2.0)*(((4*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)-(-((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)+((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)+((x2-x3)^2+(y2-y3)^2+(z2-z3)^2))^2)*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)^(-1.0)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)^(-1.0))^0.5*(2.0*(x1-x3)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)^2.0+(-x2+x3)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)*(-((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)+((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)+((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)))*((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)*(-((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)+((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)+((x2-x3)^2+(y2-y3)^2+(z2-z3)^2))+((4*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)-(-((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)+((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)+((x2-x3)^2+(y2-y3)^2+(z2-z3)^2))^2)*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)^(-1.0)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)^(-1.0))^1.5*(-x1+x2)*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)^2.0*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)^3.0)*((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)^(-0.5)*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)^(-2.0)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)^(-3.0);
        dfdy2 = ((4*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)-(-((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)+((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)+((x2-x3)^2+(y2-y3)^2+(z2-z3)^2))^2)*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)^(-1.0)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)^(-1.0))^(-2.0)*(((4*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)-(-((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)+((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)+((x2-x3)^2+(y2-y3)^2+(z2-z3)^2))^2)*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)^(-1.0)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)^(-1.0))^0.5*(2.0*(y1-y3)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)^2.0+(-y2+y3)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)*(-((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)+((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)+((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)))*((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)*(-((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)+((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)+((x2-x3)^2+(y2-y3)^2+(z2-z3)^2))+((4*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)-(-((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)+((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)+((x2-x3)^2+(y2-y3)^2+(z2-z3)^2))^2)*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)^(-1.0)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)^(-1.0))^1.5*(-y1+y2)*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)^2.0*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)^3.0)*((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)^(-0.5)*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)^(-2.0)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)^(-3.0);
        dfdz2 = ((4*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)-(-((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)+((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)+((x2-x3)^2+(y2-y3)^2+(z2-z3)^2))^2)*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)^(-1.0)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)^(-1.0))^(-2.0)*(((4*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)-(-((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)+((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)+((x2-x3)^2+(y2-y3)^2+(z2-z3)^2))^2)*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)^(-1.0)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)^(-1.0))^0.5*(2.0*(z1-z3)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)^2.0+(-z2+z3)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)*(-((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)+((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)+((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)))*((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)*(-((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)+((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)+((x2-x3)^2+(y2-y3)^2+(z2-z3)^2))+((4*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)-(-((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)+((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)+((x2-x3)^2+(y2-y3)^2+(z2-z3)^2))^2)*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)^(-1.0)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)^(-1.0))^1.5*(-z1+z2)*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)^2.0*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)^3.0)*((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)^(-0.5)*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)^(-2.0)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)^(-3.0);
        dfdx3 = ((4*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)-(-((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)+((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)+((x2-x3)^2+(y2-y3)^2+(z2-z3)^2))^2)*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)^(-1.0)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)^(-1.0))^(-1.5)*((x1-x3)*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)^2.0*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)^3.0*(-((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)+((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)+((x2-x3)^2+(y2-y3)^2+(z2-z3)^2))+(x2-x3)*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)^3.0*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)^2.0*(-((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)+((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)+((x2-x3)^2+(y2-y3)^2+(z2-z3)^2))+(-2.0*x1-2.0*x2+4.0*x3)*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)^3.0*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)^3.0)*((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)^0.5*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)^(-4.0)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)^(-4.0)*(-((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)+((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)+((x2-x3)^2+(y2-y3)^2+(z2-z3)^2));
        dfdy3 = ((4*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)-(-((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)+((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)+((x2-x3)^2+(y2-y3)^2+(z2-z3)^2))^2)*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)^(-1.0)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)^(-1.0))^(-1.5)*((y1-y3)*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)^2.0*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)^3.0*(-((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)+((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)+((x2-x3)^2+(y2-y3)^2+(z2-z3)^2))+(y2-y3)*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)^3.0*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)^2.0*(-((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)+((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)+((x2-x3)^2+(y2-y3)^2+(z2-z3)^2))+(-2.0*y1-2.0*y2+4.0*y3)*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)^3.0*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)^3.0)*((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)^0.5*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)^(-4.0)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)^(-4.0)*(-((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)+((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)+((x2-x3)^2+(y2-y3)^2+(z2-z3)^2));
        dfdz3 = ((4*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)-(-((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)+((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)+((x2-x3)^2+(y2-y3)^2+(z2-z3)^2))^2)*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)^(-1.0)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)^(-1.0))^(-1.5)*((z1-z3)*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)^2.0*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)^3.0*(-((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)+((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)+((x2-x3)^2+(y2-y3)^2+(z2-z3)^2))+(z2-z3)*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)^3.0*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)^2.0*(-((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)+((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)+((x2-x3)^2+(y2-y3)^2+(z2-z3)^2))+(-2.0*z1-2.0*z2+4.0*z3)*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)^3.0*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)^3.0)*((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)^0.5*((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)^(-4.0)*((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)^(-4.0)*(-((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)+((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)+((x2-x3)^2+(y2-y3)^2+(z2-z3)^2));

        R = [dfdx1,dfdy1,dfdz1,dfdx2,dfdy2,dfdz2,dfdx3,dfdy3,dfdz3];
    end
    Av(i,j) = R;
end
[Ai,Aj] = ndgrid((1:nTrig),(1:nC));
Re = sparse(Ai(:),Aj(:),Av(:),nTrig,nC)*dNd;
end







 