function [KeL] = KeL_Frame(Coor,SecProp,MatProp,varargin)
A = SecProp(1); J = SecProp(2); Iy = SecProp(3); Iz = SecProp(4); 
E = MatProp(1); G = MatProp(2);
Xt = (Coor(2,:)-Coor(1,:)); 
 L  = norm(Xt); 
if nargin < 4
    Linv = @(x) L^-x;
elseif strcmp(varargin{1},'shape')
    d   = varargin{2};
    dCs = d==reshape((1:9),3,3);
    dXt = (dCs(2,:)-dCs(1,:));
    dL  = dot(Xt,dXt)/L;
    Linv = @(x) -x*L^-(x+1)*dL;
else 
    Linv = @(x) L^-x;
end
M1 = [1,-1;-1,1];
M2 = [0, 1;-1,0];
KeL = zeros(12);
    KeL([1       7],[1       7]) = E*A*Linv(1)*M1;
    KeL([4      10],[4      10]) = G*J*Linv(1)*M1;
    KeL([2 3  8  9],[2 3  8  9]) = 12*E*Linv(3)*[eye(2),-eye(2);-eye(2),eye(2)];
    KeL([5 6 11 12],[5 6 11 12]) = 2*E*Linv(1)*[2*eye(2),eye(2);eye(2),2*eye(2)];
    KeL([2 3  8  9],[5 6 11 12]) = 6*E*Linv(2)*[M2,M2;-M2,-M2];
    KeL([5 6 11 12],[2 3  8  9]) = 6*E*Linv(2)*[-M2,M2;-M2,M2];
    KeL([3 5  9 11],:          ) = Iy*KeL([3 5 9 11],:);
    KeL([2 6  8 12],:          ) = Iz*KeL([2 6 8 12],:);
end

% M1 = [1,-1;-1,1];
% M2 = [2  1; 1 2];
% KeL = zeros(12);
%     KeL([1       7],[1       7]) = E*A*Linv(1)*M1;
%     KeL([4      10],[4      10]) = G*J*Linv(1)*M1;
%     KeL([2 3  8  9],[2 3  8  9]) = 12*E*Linv(3)*[M1,0*M1;0*M1,M1];
%     KeL([5 6 11 12],[5 6 11 12]) =  2*E*Linv(1)*[M2,0*M2;0*M2,M2];
%     KeL([2 3  8  9],[5 6 11 12]) = kron([1 1;-1 -1],6*E*Linv(2)*[0 1;-1 0]);
%     KeL([5 6 11 12],[2 3  8  9]) = kron([1 -1;1 -1],6*E*Linv(2)*[0 -1;1 0]);
%     KeL([3 5  9 11],:          ) = Iy*KeL([3 5 9 11],:);
%     KeL([2 6  8 12],:          ) = Iz*KeL([2 6 8 12],:);






% KeL = zeros(12);
%     KeL([1       7],[1       7]) = E*A*Linv(1)*[1 -1; -1 1];
%     KeL([4      10],[4      10]) = G*J*Linv(1)*[1 -1; -1 1];
%     KeL([2 3  8  9],[2 3  8  9]) = kron(12*E*Linv(3)*[1 -1;-1 1],eye(2));
%     KeL([5 6 11 12],[5 6 11 12]) = kron( 2*E*Linv(1)*[2  1; 1 2],eye(2));
%     KeL([2 3  8  9],[5 6 11 12]) = kron([1 1;-1 -1],6*E*Linv(2)*[0 1;-1 0]);
%     KeL([5 6 11 12],[2 3  8  9]) = kron([1 -1;1 -1],6*E*Linv(2)*[0 -1;1 0]);
%     KeL([3 5  9 11],:          ) = Iy*KeL([3 5 9 11],:);
%     KeL([2 6  8 12],:          ) = Iz*KeL([2 6 8 12],:);