function [Ke] = KeLCS_TrigMembrane(locC, h, E, nu)
% This function generates the element stiffness matrix of a membrane 
% element wrt a local coordinate system. The DOFs are arranged per 
% node as [u, v], where u is the displacement wrt the local x-axis 
% and v is the displacement wrt the local y-axis.
%    Input:   locC = the x- and y-coordinates of each vertex wrt a
%                    local coordinate system [x1 y1; x2 y2; x3 y3]    
%             h    = the thickness of the plate
%             E    = the Young’s modulus of the plate’s material
%             nu   = the Poisson’s ratio of the plate’s material
%    Output:  Ke   = The locally defined element stiffness matrix

Shape = trigShape(locC);
Area  = trigSArea(locC);
t     = sum(h(:),1)/size(h(:),1);
D     = E*[1,nu,0;nu,1,0;0,0,(1-nu)/2]/(1-nu^2);
B     = kron(Shape(:,1)',[1,0;0,0;0,1])+...
        kron(Shape(:,2)',[0,0;0,1;1,0]);
Ke    = B'*D*B*Area*t;
end

% input:  locally defined coordinates [x1 y1; x2 y2; x3 y3]
% output: shape functions [N1; N2; N3] with Ni = Xc,i*x + Yc,i*y + Cc,i
function [Shape] = trigShape(locC)
X  = locC(:,1);                 Y  = locC(:,2);
Xc = Y([2 3 1])-Y([3 1 2]);     Yc = X([3 1 2])-X([2 3 1]);
Cc = X([2 3 1]).*Y([3 1 2])-X([3 1 2]).*Y([2 3 1]);
AA = (X(2)-X(1))*(Y(3)-Y(1))-(X(3)-X(1))*(Y(2)-Y(1));
Shape = [Xc Yc Cc]/AA;
end
% input:  locally defined coordinates [x1 x2 x3; y1 y2 y3]
% output: surface area of triangular element
function [SurfArea] = trigSArea(locC)
X  = locC(:,1);                 Y  = locC(:,2);
SurfArea = 0.5*det([X(2:3)-X(1), Y(2:3)-Y(1)]);
end