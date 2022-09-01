function [out]=C_Circle_Trig(CL,varargin)
Solv = @(x,y) det3([x,y,ones(3,1)]);
CL2 = sum(CL.^2,2);
Nom = Solv(CL(:,1),CL(:,2));
Nx0  = Solv(CL2,CL(:,2));
Ny0  = Solv(CL(:,1),CL2);
if nargin == 1
    x0 = Nx0/Nom/2; y0 = Ny0/Nom/2;
elseif nargin > 1
   dCL  = varargin{2}==reshape((1:9),3,3);
   dCL2 = sum(2*CL.*dCL,2);
   dFrac = @(t,dt,n,dn) (dt*n-t*dn)*n^-2;
   dNom = Solv(dCL(:,1),CL(:,2))+Solv(CL(:,1),dCL(:,2));
   dNx0  = Solv(dCL2,CL(:,2))+Solv(CL2,dCL(:,2));
   dNy0  = Solv(dCL(:,1),CL2)+Solv(CL(:,1),dCL2);
   x0 = dFrac(Nx0,dNx0,Nom,dNom)/2;
   y0 = dFrac(Ny0,dNy0,Nom,dNom)/2;
end
out = [x0,y0,0];
end

