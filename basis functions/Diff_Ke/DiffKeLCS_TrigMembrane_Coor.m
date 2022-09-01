
function [dKdx1,dKdy1,dKdx2,dKdy2,dKdx3,dKdy3] = ...
          DiffKeLCS_TrigMembrane_Coor(locC, h, E, nu)
% This function differentiates the element stiffness matrix of a membrane 
% element wrt the locally defined vertex coordinates. The DOFs are arranged 
% per node as [u, v], where u is the displacement wrt the local x-axis and 
% v is the displacement wrt the local y-axis.
%    Input:   locC  = the x- and y-coordinates of each vertex wrt a
%                     local coordinate system [x1 y1; x2 y2; x3 y3]    
%             h     = the thickness of the plate
%             E     = the Young’s modulus of the plate’s material
%             nu    = the Poisson’s ratio of the plate’s material
%    Output:  dKdxi = A 6x6 matrix containing the derivatives of each 
%                     element of the locally defined element stiffness 
%                     matrix wrt the local x-coordinate of vertex i
%             dKdyi = A 6x6 matrix containing the derivatives of each 
%                     element of the locally defined element stiffness 
%                     matrix wrt the local y-coordinate of vertex i

t     = sum(h(:),1)/size(h(:),1);
D     = E*[1,nu,0;nu,1,0;0,0,(1-nu)/2]/(1-nu^2);
dKdx1 = diffKe_locC(locC,D,t,1);  dKdy1 = diffKe_locC(locC,D,t,4);
dKdx2 = diffKe_locC(locC,D,t,2);  dKdy2 = diffKe_locC(locC,D,t,5);
dKdx3 = diffKe_locC(locC,D,t,3);  dKdy3 = diffKe_locC(locC,D,t,6);
end

function [dK] = diffKe_locC(locC,D,t,d)
[Shape,dShape,Area,dArea] = diffTrigShapeArea_Coor(locC,d);
H1    = [1,0;0,0;0,1];    H2    = [0,0;0,1;1,0];
B     = kron( Shape(:,1)',H1)+kron( Shape(:,2)',H2);
dB    = kron(dShape(:,1)',H1)+kron(dShape(:,2)',H2);
dK    = (dB'*D*B*Area+B'*D*dB*Area+B'*D*B*dArea)*t;
end

% input:  locC: locally defined coordinates [x1 y1; x2 y2; x3 y3]
%         d:    index of differentiation variable e.g. x2=2, y3=6
% output: Shape functions, area and their derivatives wrt given coordinate
function [Shape,dShape,Area,dArea] = diffTrigShapeArea_Coor(locC,d)
X   = locC(:,1);                 dX   = [d==1;d==2;d==3];          
Y   = locC(:,2);                 dY   = [d==4;d==5;d==6];
Xc  =  Y([2 3 1])- Y([3 1 2]);   dXc  = dY([2 3 1])-dY([3 1 2]);
Yc  =  X([3 1 2])- X([2 3 1]);   dYc  = dX([3 1 2])-dX([2 3 1]);
AA  = det([ X(2:3)- X(1), Y(2:3)- Y(1)]);
dAA = det([dX(2:3)-dX(1), Y(2:3)- Y(1)])+...
      det([ X(2:3)- X(1),dY(2:3)-dY(1)]);
Shape = [Xc Yc]/AA;      dShape = ([dXc dYc]*AA-[Xc Yc]*dAA)*AA^-2;
Area  = 0.5*AA;          dArea  = 0.5*dAA;
end