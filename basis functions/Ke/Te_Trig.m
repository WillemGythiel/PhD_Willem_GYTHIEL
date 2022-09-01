function [out] = Te_Trig(Cor,varargin)
   D =  Cor-repmat( Cor(1,:),3,1);
unit = @(x) {x/norm(x),x,norm(x)};
   Z = unit(cross3(D(2,:),D(3,:)));
   X = unit(D(2,:));
   Y = unit(cross3(Z{1},X{1}));
  Te = [X{1};Y{1};Z{1}]; 
if nargin == 1
 out = Te;
elseif nargin > 1
   dCor = varargin{2}==reshape((1:9),3,3);
     dD = dCor-repmat(dCor(1,:),3,1);
   diff = @(v,a) (a/v{3}-v{2}*a'*v{2}/v{3}^3);  
    deZ = diff(Z,cross3(dD(2,:),D(3,:))+cross3(D(2,:),dD(3,:)));
    deX = diff(X,dD(2,:));
    deY = diff(Y,cross3(deZ,X{1})+cross3(Z{1},deX));
    out = [deX;deY;deZ]; 
end
end
