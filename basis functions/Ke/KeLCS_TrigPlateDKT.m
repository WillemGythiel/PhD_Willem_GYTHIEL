function [Ke] = KeLCS_TrigPlateDKT(locC, h, E, nu)
% This function generates the element stiffness matrix of a DKT plate 
% element wrt a local coordinate system. The DOFs are arranged per node 
% as [w, theta_x and theta_y], where w is the deflection, theta_x is the 
% derivative of the deflection wrt x (hence oriented along the y-axis) 
% and theta_y is the derivative of the deflection wrt y (hence oriented 
% along the negative x-axis), as is conventional for in plane rotations. 
%    Input:   locC = the x- and y-coordinates of each vertex wrt a
%                    local coordinate system [x1 y1; x2 y2; x3 y3]    
%             h    = the thickness of the plate
%             E    = the Young’s modulus of the plate’s material
%             nu   = the Poisson’s ratio of the plate’s material
%    Output:  Ke   = The locally defined element stiffness matrix

Da    = E*h^3*[1,nu,0;nu,1,0;0,0,(1-nu)/2]/(12*(1-nu^2));
Shape = trigShape(locC);
AA    = trigSAA(locC);
Ka    = KpMaker(Da,Shape,AA);
Kb    = KpMaker(Da,Shape+Shape([2 3 1],:),AA);
CofM  = kron(-(1/6)*ones(3)+(2/3)*eye(3),ones(2));
  K1  = CofM.*Ka;
  K2  = 2/3*K2Maker(Ka);
  K3  = 2/3*(Kb + AuxMaker(Ka));
Krot  = [K1 K2'; K2 K3];
Ts    = TsMat(locC);
Ke    = Ts'*Krot*Ts;  
end
% input:  locally defined coordinates [x1 y1; x2 y2; x3 y3]
% output: shape functions [N1; N2; N3] with Ni = Xc,i*x + Yc,i*y + Cc,i
function [Shape]  = trigShape(locC)
X  = locC(:,1);                 Y  = locC(:,2);
Xc = Y([2 3 1])-Y([3 1 2]);     Yc = X([3 1 2])-X([2 3 1]);
Cc = X([2 3 1]).*Y([3 1 2])-X([3 1 2]).*Y([2 3 1]);
AA = (X(2)-X(1))*(Y(3)-Y(1))-(X(3)-X(1))*(Y(2)-Y(1));
Shape = [Xc Yc Cc]/AA;
end
% input:  locally defined coordinates [x1 y1; x2 y2; x3 y3]
% output: 2*surface area of triangular element
function [Area2]  = trigSAA(locC)
X  = locC(:,1);                 Y  = locC(:,2);
Area2 = det([X(2:3)-X(1), Y(2:3)-Y(1)]);
end
function [K]      = KpMaker(Da,Shape,AA)
H1    = [1,0;0,0;0,1];    H2    = [0,0;0,1;1,0];
B     = kron( Shape(:,1)',H1)+kron( Shape(:,2)',H2);
K     = B'*Da*B*AA;
end
function [K2]     = K2Maker(Kinput)
K2 = 0*Kinput;
  K2(1:2,3:4) = Kinput(1:2,3:4);
  K2(3:4,5:6) = Kinput(3:4,5:6);
  K2(5:6,1:2) = Kinput(5:6,1:2);
  K2(1:2,1:2) = Kinput(3:4,1:2);
  K2(3:4,3:4) = Kinput(5:6,3:4);
  K2(5:6,5:6) = Kinput(1:2,5:6);
end
function [Auxmat] = AuxMaker(Kinput)
Auxmat = 0*Kinput;
  Auxmat(1:2,1:2) = Kinput(1:2,1:2)+Kinput(3:4,3:4);
  Auxmat(3:4,3:4) = Kinput(3:4,3:4)+Kinput(5:6,5:6);
  Auxmat(5:6,5:6) = Kinput(5:6,5:6)+Kinput(1:2,1:2);
  Auxmat(3:4,1:2) = Kinput(5:6,1:2);
  Auxmat(5:6,1:2) = Kinput(5:6,3:4);
  Auxmat(1:2,3:4) = Kinput(1:2,5:6);
  Auxmat(5:6,3:4) = Kinput(1:2,3:4);
  Auxmat(1:2,5:6) = Kinput(3:4,5:6);
  Auxmat(3:4,5:6) = Kinput(3:4,1:2);
end
function [TsM]    = TsMat(locC)
TsM                    = zeros(12,9);
TsM(1:6,[2:3,5:6,8:9]) = eye(6); 
TsM( 7: 8, 1:6 )       = CoefTs(locC([1 2],:));
TsM( 9:10, 4:9 )       = CoefTs(locC([2 3],:));
TsM(11:12,[7:9,1:3])   = CoefTs(locC([3 1],:));
end
function [Res]    = CoefTs(Cs)
Xij = Cs(2)-Cs(1); Yij =  Cs(4)- Cs(3);         Num = 4*(Xij^2+Yij^2);
% the coefficients
Aij = 6*Xij/Num;   Bij = (2*Yij^2-Xij^2)/Num;   Cij = -3*Xij*Yij/Num;
Dij = 6*Yij/Num;   Eij = (2*Xij^2-Yij^2)/Num;
% the output
Res  = [-Aij,Bij,Cij,Aij,Bij,Cij;...
        -Dij,Cij,Eij,Dij,Cij,Eij];
end