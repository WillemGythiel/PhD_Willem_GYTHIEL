function [dKdx1,dKdy1,dKdx2,dKdy2,dKdx3,dKdy3] = ...
          DiffKeLCS_TrigPlateDKT_Coor(locC, h, E, nu)
% This function differentiates the element stiffness matrix of a DKT 
% plate element wrt the locally defined vertex coordinates. The DOFs 
% are arranged per node as [w, theta_x and theta_y], where w is the 
% deflection, theta_x is the derivative of the deflection wrt x (hence 
% oriented along the y-axis) and theta_y is the derivative of the 
% deflection wrt y (hence oriented along the negative x-axis), as is 
% conventional for in plane rotations. 
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

Da    = E*h^3*[1,nu,0;nu,1,0;0,0,(1-nu)/2]/(12*(1-nu^2));
dKdx1 = diffKe_locC(locC,Da,1);  dKdy1 = diffKe_locC(locC,Da,4);
dKdx2 = diffKe_locC(locC,Da,2);  dKdy2 = diffKe_locC(locC,Da,5);
dKdx3 = diffKe_locC(locC,Da,3);  dKdy3 = diffKe_locC(locC,Da,6);
end
function [dKe]    = diffKe_locC(locC,Da,d)
Ts       = TsMat(locC);
dTs      = diffTsMat_Coor(locC,d);
[Kr,dKr] = diffKrot_Coor(locC,Da,d);
dKe      = dTs'*Kr*Ts+Ts'*dKr*Ts+Ts'*Kr*dTs;
end
function [Kr,dKr] = diffKrot_Coor(locC,Da,d)
[Shape,dShape,AA,dAA] = diffTrigShapeAA_Coor(locC,d);
[Ka,dKa] = dKpMaker(Da, Shape,dShape,AA,dAA);
[Kb,dKb] = dKpMaker(Da, Shape+ Shape([2 3 1],:),...
                       dShape+dShape([2 3 1],:),AA,dAA);
Kr       = KrotMaker(Ka,Kb);
dKr      = KrotMaker(dKa,dKb);
end
function [Krot]   = KrotMaker(Ka,Kb)
CofM = kron(-(1/6)*ones(3)+(2/3)*eye(3),ones(2));
  K1 = CofM.*Ka;
  K2 = 2/3*K2Maker(Ka);
  K3 = 2/3*(Kb + AuxMaker(Ka));
Krot = [K1 K2'; K2 K3];
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
function [K,dK]   = dKpMaker(Da,Shape,dShape,AA,dAA)
H1    = [1,0;0,0;0,1];    H2    = [0,0;0,1;1,0];
B     = kron( Shape(:,1)',H1)+kron( Shape(:,2)',H2);
dB    = kron(dShape(:,1)',H1)+kron(dShape(:,2)',H2);
K     = B'*Da*B*AA;
dK    = dB'*Da*B*AA+B'*Da*dB*AA+B'*Da*B*dAA;
end
function [Shape,dShape,AA,dAA] = diffTrigShapeAA_Coor(locC,d)
X   = locC(:,1);                 dX   = [d==1;d==2;d==3];          
Y   = locC(:,2);                 dY   = [d==4;d==5;d==6];
Xc  =  Y([2 3 1])- Y([3 1 2]);   dXc  = dY([2 3 1])-dY([3 1 2]);
Yc  =  X([3 1 2])- X([2 3 1]);   dYc  = dX([3 1 2])-dX([2 3 1]);
AA  = det([ X(2:3)- X(1), Y(2:3)- Y(1)]);
dAA = det([dX(2:3)-dX(1), Y(2:3)- Y(1)])+...
      det([ X(2:3)- X(1),dY(2:3)-dY(1)]);
Shape = [Xc Yc]/AA;      dShape = ([dXc dYc]*AA-[Xc Yc]*dAA)*AA^-2;
end
function [TsM] = TsMat(locC)
TsM                    = zeros(12,9);
TsM(1:6,[2:3,5:6,8:9]) = eye(6); 
TsM( 7: 8, 1:6 )       = CoefTs(locC([1 2],:));
TsM( 9:10, 4:9 )       = CoefTs(locC([2 3],:));
TsM(11:12,[7:9,1:3])   = CoefTs(locC([3 1],:));
end
function [Res] = CoefTs(Cs)
Xij = Cs(2)-Cs(1); Yij =  Cs(4)- Cs(3);         Num = 4*(Xij^2+Yij^2);
% the coefficients
Aij = 6*Xij/Num;   Bij = (2*Yij^2-Xij^2)/Num;   Cij = -3*Xij*Yij/Num;
Dij = 6*Yij/Num;   Eij = (2*Xij^2-Yij^2)/Num;
% the output
Res  = [-Aij,Bij,Cij,Aij,Bij,Cij;...
        -Dij,Cij,Eij,Dij,Cij,Eij];
end
function [dTs] = diffTsMat_Coor(locC,d)
dLocC                = [d==1,d==4;d==2,d==5;d==3,d==6];
dTs                  = zeros(12,9);
dTs( 7: 8, 1:6 )     = diffCoefTs_Coor(locC([1 2],:),dLocC([1 2],:));
dTs( 9:10, 4:9 )     = diffCoefTs_Coor(locC([2 3],:),dLocC([2 3],:));
dTs(11:12,[7:9,1:3]) = diffCoefTs_Coor(locC([3 1],:),dLocC([3 1],:));
end
function [Res] = diffCoefTs_Coor(Cor,dCor)
Xij  =  Cor(2)- Cor(1);      dXij  = dCor(2)-dCor(1);   
Yij  =  Cor(4)- Cor(3);      dYij  = dCor(4)-dCor(3);
Num  = 4*(Xij^2+Yij^2);      dNum  = 8*(Xij*dXij+Yij*dYij);
% the derivatives of the coefficients
dAij = (( 6*dXij                )*Num-( 6*Xij        )*dNum)*Num^-2;
dBij = (( 4*Yij*dYij-2*Xij*dXij )*Num-( 2*Yij^2-Xij^2)*dNum)*Num^-2;
dCij = ((-3*(dXij*Yij+Xij*dYij) )*Num-(-3*Xij*Yij    )*dNum)*Num^-2;
dDij = (( 6*dYij                )*Num-( 6*Yij        )*dNum)*Num^-2;
dEij = (( 4*Xij*dXij-2*Yij*dYij )*Num-( 2*Xij^2-Yij^2)*dNum)*Num^-2;
% the output
Res  = [-dAij,dBij,dCij,dAij,dBij,dCij;...
        -dDij,dCij,dEij,dDij,dCij,dEij];
end