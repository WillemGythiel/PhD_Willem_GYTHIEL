function [out] = Coords_Trig(Cor,varargin)
   D  =  Cor-repmat( Cor(1,:),3,1);
   Te = Te_Trig(Cor)';
if nargin == 1
   out = D*Te;
elseif nargin > 1
   dCor = varargin{2}==reshape((1:9),3,3);
   dD   = dCor-repmat(dCor(1,:),3,1);
   dTe  = Te_Trig(Cor,'shape',varargin{2})';
   out  = dD*Te+D*dTe;
end  
end