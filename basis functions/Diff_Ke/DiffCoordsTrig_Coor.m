function [dCLdx1,dCLdy1,dCLdz1,dCLdx2,dCLdy2,dCLdz2,dCLdx3,dCLdy3,dCLdz3]=...
          DiffCoordsTrig_Coor(Cs)
dCLdx1 = diffCL(Cs,1);  dCLdy1 = diffCL(Cs,4);  dCLdz1 = diffCL(Cs,7);
dCLdx2 = diffCL(Cs,2);  dCLdy2 = diffCL(Cs,5);  dCLdz2 = diffCL(Cs,8);
dCLdx3 = diffCL(Cs,3);  dCLdy3 = diffCL(Cs,6);  dCLdz3 = diffCL(Cs,9);
end
function [dCdsLoc]=diffCL(Coords,d)
[Te,dTe] = diffTeMatrix(Coords,d);
dCoords  = d==reshape((1:9),3,3);
      D  =  Coords-repmat( Coords(1,:),3,1);
     dD  = dCoords-repmat(dCoords(1,:),3,1);
dCdsLoc  = dD*Te'+D*dTe';
end
function [Te,dTe]=diffTeMatrix(Coords,d)
dCoords = d==reshape((1:9),3,3);
%-----------------------------------------------------------------------e_z
 NumerZ  = cross( Coords(2,:)- Coords(1,:), Coords(3,:)- Coords(1,:));
dNumerZ  = cross(dCoords(2,:)-dCoords(1,:), Coords(3,:)- Coords(1,:))+...
           cross( Coords(2,:)- Coords(1,:),dCoords(3,:)-dCoords(1,:));
[eZ,deZ] = diffUnitVector(NumerZ, dNumerZ);
%-----------------------------------------------------------------------e_x
 NumerX  = 0.5*( Coords(2,:)+ Coords(3,:))- Coords(1,:);
dNumerX  = 0.5*(dCoords(2,:)+dCoords(3,:))-dCoords(1,:);
[eX,deX] = diffUnitVector(NumerX, dNumerX);
%-----------------------------------------------------------------------e_y
 NumerY  = cross(eZ,eX);
dNumerY  = cross(deZ,eX)+cross(eZ,deX);
[eY,deY] = diffUnitVector(NumerY, dNumerY);
Te = [ eX; eY; eZ];    dTe = [deX;deY;deZ];
end
function [UnitV,dUnitV]=diffUnitVector(Numer, dNumer)
Denom = norm(Numer);    dDenom =  Numer*dNumer'/Denom;
UnitV = Numer/Denom;    dUnitV = (dNumer*Denom-Numer*dDenom)*Denom^-2;
end


% TeZx = (y2-y1)*(z3-z1)-(y3-y1)*(z2-z1);
% TeZy = (x3-x1)*(z2-z1)-(x2-x1)*(z3-z1);
% TeZz = (x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);
% TeZ  = [TeZx TeZy TeZz];
% % denominator of e_z = 2Delta
% NeZ = (TeZx^2 + TeZy^2 + TeZz^2)^0.5;
% % complete formulation of e_z
% eZ = TeZ/NeZ;
% % derivative of the numerator of e_z
% dTeZx = (dy2-dy1)*( z3- z1)-(dy3-dy1)*( z2- z1)+...
%         ( y2- y1)*(dz3-dz1)-( y3- y1)*(dz2-dz1);
% dTeZy = (dx3-dx1)*( z2- z1)-(dx2-dx1)*( z3- z1)+...
%         ( x3- x1)*(dz2-dz1)-( x2- x1)*(dz3-dz1);
% dTeZz = (dx2-dx1)*( y3- y1)-(dx3-dx1)*( y2- y1)+...
%         ( x2- x1)*(dy3-dy1)-( x3- x1)*(dy2-dy1);
% dTeZ = [dTeZx dTeZy dTeZz];