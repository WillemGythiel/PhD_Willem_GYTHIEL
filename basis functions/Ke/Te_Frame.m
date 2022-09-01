function [Te] = Te_Frame(Coor,varargin)
Vt = (Coor(3,:)-Coor(1,:));
Xt = (Coor(2,:)-Coor(1,:));
Yt = cross3(Vt,Xt);
Zt = cross3(Xt,Yt);
if nargin < 2
    Tm  = [Xt/norm(Xt);Yt/norm(Yt);Zt/norm(Zt)];
elseif strcmp(varargin{1},'size')
    Tm  = zeros(3);
elseif strcmp(varargin{1},'shape')
    d   = varargin{2};
    dCs = d==reshape((1:9),3,3);
    dVt = (dCs(3,:)-dCs(1,:));
    dXt = (dCs(2,:)-dCs(1,:));
    dYt = cross3(dVt,Xt)+cross3(Vt,dXt);
    dZt = cross3(dXt,Yt)+cross3(Xt,dYt);
    Tm  = [dUV(Xt,dXt);dUV(Yt,dYt);dUV(Zt,dZt)];
end
Te  = kron(eye(4),Tm);
end
function [dUnitV]=dUV(Numer,dNumer)
Denom = norm(Numer);    dDenom =  dot3(Numer,dNumer)/Denom;
dUnitV = (dNumer*Denom-Numer*dDenom)*Denom^-2;
end