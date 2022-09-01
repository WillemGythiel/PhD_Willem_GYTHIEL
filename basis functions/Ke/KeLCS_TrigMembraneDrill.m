function [Ke] = KeLCS_TrigMembraneDrill(locC, h, E, nu)
% This function generates the element stiffness matrix of a membrane 
% element wrt a local coordinate system. The DOFs are arranged per 
% node as [u, v, theta_z], where u is the displacement wrt the local 
% x-axis, v is the displacement wrt the local y-axis and theta_z is 
% the rotation wrt the local z-axis (the drilling DOF).
%    Input:   locC = the x- and y-coordinates of each vertex wrt a
%                    local coordinate system [x1 y1; x2 y2; x3 y3]    
%             h    = the thickness of the plate
%             E    = the Young’s modulus of the plate’s material
%             nu   = the Poisson’s ratio of the plate’s material
%    Output:  Ke   = The locally defined element stiffness matrix

Shape = trigShape(locC);        Area  = trigSArea(locC);
X     = locC(:,1);              Y     = locC(:,2);
xT    = sum(X)/3-X';            yT    = sum(Y)/3-Y';
D     = E*[1,nu,0;nu,1,0;0,0,(1-nu)/2]/(1-nu^2);
H1    = [1,0,0;0,0,0;0,1,0];    H2    = [0,0,0;0,1,0;1,0,0];
B     = kron(Shape(:,1)',H1)+kron(Shape(:,2)',H2)+...
        kron([[-yT;xT].*Shape';sum([xT;-yT].*Shape',1)],[0,0,1]);
K1    = KauxMaker(sum(X),sum(Y),Shape,D)/-36;
K2    = KauxMaker(X(1),Y(1),Shape,D)/12;
K3    = KauxMaker(X(2),Y(2),Shape,D)/12; 
K4    = KauxMaker(X(3),Y(3),Shape,D)/12; 
Ke    = (B'*D*B+K1+K2+K3+K4)*Area*h;
end

% input:  locally defined coordinates [x1 y1; x2 y2; x3 y3]
% output: shape functions [N1; N2; N3] with Ni = Xc,i*x + Yc,i*y + Cc,i
function [Shape] = trigShape(locC)
X  = locC(:,1);                 Y  = locC(:,2);
Xc = Y([2 3 1])-Y([3 1 2]);     Yc = X([3 1 2])-X([2 3 1]);
% Cc = X([2 3 1]).*Y([3 1 2])-X([3 1 2]).*Y([2 3 1]);
AA = (X(2)-X(1))*(Y(3)-Y(1))-(X(3)-X(1))*(Y(2)-Y(1));
Shape = [Xc Yc]/AA;
end
% input:  locally defined coordinates [x1 y1; x2 y2; x3 y3]
% output: surface area of triangular element
function [SurfArea] = trigSArea(locC)
X  = locC(:,1);                 Y  = locC(:,2);
SurfArea = 0.5*det([X(2:3)-X(1), Y(2:3)-Y(1)]);
end

function [Kc] = KauxMaker(Xi,Yi,Shape,D)
Bc = kron([[-Yi; Xi ].*Shape';sum([ Xi;-Yi ].*Shape',1)],[0,0,1]);
Kc = Bc'*D*Bc;
end