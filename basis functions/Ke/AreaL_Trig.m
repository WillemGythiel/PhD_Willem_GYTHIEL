function [out]=AreaL_Trig(Cor,varargin)
   D = Cor-repmat(Cor(1,:),3,1);
 Num = cross3(D(2,:),D(3,:));
if nargin == 1
 out = Num(3)/2;
elseif nargin > 1
   dCor = varargin{2}==reshape((1:9),3,3);
     dD = dCor-repmat(dCor(1,:),3,1);
   dNum = cross3(dD(2,:),D(3,:))+cross3(D(2,:),dD(3,:));
    out = dNum(3)/2;
end
end