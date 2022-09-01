
function [dKdx1,dKdy1,dKdx2,dKdy2,dKdx3,dKdy3] = ...
          DiffKeLCS_TrigMembraneDrill_Coor(locC, h, E, nu)
% This function differentiates the element stiffness matrix of a membrane 
% element wrt the locally defined vertex coordinates. The DOFs are arranged 
% per node as [u, v, theta_z], where u is the displacement wrt the local 
% x-axis, v is the displacement wrt the local y-axis and theta_z is the 
% rotation wrt the local z-axis (the drilling DOF).
%    Input:   locC  = the x- and y-coordinates of each vertex wrt a
%                     local coordinate system [x1 y1; x2 y2; x3 y3]    
%             h     = the thickness of the plate
%             E     = the Young’s modulus of the plate’s material
%             nu    = the Poisson’s ratio of the plate’s material
%    Output:  dKdxi = A 9x9 matrix containing the derivatives of each 
%                     element of the locally defined element stiffness 
%                     matrix wrt the local x-coordinate of vertex i
%             dKdyi = A 9x9 matrix containing the derivatives of each 
%                     element of the locally defined element stiffness 
%                     matrix wrt the local y-coordinate of vertex i

D     = E*[1,nu,0;nu,1,0;0,0,(1-nu)/2]/(1-nu^2);
dKdx1 = diffKe_locC(locC,D,h,1);  dKdy1 = diffKe_locC(locC,D,h,4);
dKdx2 = diffKe_locC(locC,D,h,2);  dKdy2 = diffKe_locC(locC,D,h,5);
dKdx3 = diffKe_locC(locC,D,h,3);  dKdy3 = diffKe_locC(locC,D,h,6);
end

function [dK] = diffKe_locC(locC,D,t,d)
X   = locC(:,1);              dX   = [d==1;d==2;d==3];          
Y   = locC(:,2);              dY   = [d==4;d==5;d==6];
[Shp,dShp,Area,dArea] = diffTrigShapeArea_Coor(X,dX,Y,dY);
[Kp, dKp] = KMaker(X,dX,Y,dY,Shp,dShp,D);
[Kc, dKc] = KauxMaker(X,dX,Y,dY,Shp,dShp,D);
dK  = ((dKp+dKc)*Area+(Kp+Kc)*dArea)*t;
end

% input:  locC: locally defined coordinates [x1 y1; x2 y2; x3 y3]
%         d:    index of differentiation variable e.g. x2=2, y3=6
% output: Shape functions, area and their derivatives wrt given coordinate
function [Shape,dShape,Area,dArea] = diffTrigShapeArea_Coor(X,dX,Y,dY)
Xc  =  Y([2 3 1])- Y([3 1 2]);   dXc  = dY([2 3 1])-dY([3 1 2]);
Yc  =  X([3 1 2])- X([2 3 1]);   dYc  = dX([3 1 2])-dX([2 3 1]);
AA  = det([ X(2:3)- X(1), Y(2:3)- Y(1)]);
dAA = det([dX(2:3)-dX(1), Y(2:3)- Y(1)])+...
      det([ X(2:3)- X(1),dY(2:3)-dY(1)]);
Shape = [Xc Yc]/AA;      dShape = ([dXc dYc]*AA-[Xc Yc]*dAA)*AA^-2;
Area  = 0.5*AA;          dArea  = 0.5*dAA;
end

function [NxT] = NxTil(xT,yT,Shape)
NxT  = [[-yT;xT].*Shape';sum([xT;-yT].*Shape',1)];
end
function [Kp, dKp] = KMaker(X,dX,Y,dY,Shp,dShp,D)
xT  = sum(X)/3-X';            dxT  = sum(dX)/3-dX';
yT  = sum(Y)/3-Y';            dyT  = sum(dY)/3-dY';
NxT = NxTil( xT, yT,Shp);     dNxT = NxTil(dxT,dyT,Shp)+NxTil(xT,yT,dShp);
H1  = [1,0,0;0,0,0;0,1,0];    H2   = [0,0,0;0,1,0;1,0,0]; H3 = [0,0,1];
B   = kron( Shp(:,1)',H1)+kron( Shp(:,2)',H2)+kron( NxT,H3);
dB  = kron(dShp(:,1)',H1)+kron(dShp(:,2)',H2)+kron(dNxT,H3);
Kp  = B'*D*B;                 dKp  = dB'*D*B+B'*D*dB;
end
function [Kaux, dKaux] = KauxMaker(X,dX,Y,dY,Shp,dShp,D)
B1 = Baux(sum(X),sum(Y),Shp);  dB1 = (Baux(sum(dX),sum(dY),Shp)+Baux(sum(X),sum(Y),dShp));  
B2 = Baux(X(1),Y(1),Shp);      dB2 = (Baux(dX(1),dY(1),Shp)+Baux(X(1),Y(1),dShp));
B3 = Baux(X(2),Y(2),Shp);      dB3 = (Baux(dX(2),dY(2),Shp)+Baux(X(2),Y(2),dShp));      
B4 = Baux(X(3),Y(3),Shp);      dB4 = (Baux(dX(3),dY(3),Shp)+Baux(X(3),Y(3),dShp)); 
Kaux  = ( B1'*D*B1/-3 + B2'*D*B2     +  B3'*D*B3 + B4'*D*B4)/12;
dKaux = (dB1'*D*B1/-3 + B1'*D*dB1/-3 + dB2'*D*B2 + B2'*D*dB2+...
         dB3'*D*B3    + B3'*D*dB3    + dB4'*D*B4 + B4'*D*dB4)/12;
end
function [Bc] = Baux(Xi,Yi,Shape)
Bc = kron([[-Yi; Xi ].*Shape';sum([ Xi;-Yi ].*Shape',1)],[0,0,1]);
end


