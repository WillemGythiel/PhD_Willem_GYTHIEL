function [K] = KtrussMaker(Nds,Brs,E,A)
nDOF = 2*size(Nds,1); 
nBrI = size(Brs,1);
Mi = zeros(4*nBrI,1); Mj = Mi;  Mv = Mi;
for i = 1:nBrI
    C21 = Nds(Brs(i,2),:)-Nds(Brs(i,1),:);
    Lt = norm(C21);
    Tm = kron(eye(2),C21(:,1:2))/Lt;
    Kl = E*A(i)/Lt*[1,-1;-1,1];
    Ke = Tm'*Kl*Tm;
    [Ai,Aj] = ndgrid(Addkron(2*(Brs(i,:)-1),(1:2)));
    ind = 16*(i-1)+(1:16);
    Mi(ind) = Ai(:);
    Mj(ind) = Aj(:);
    Mv(ind) = Ke(:);
end
K = sparse(Mi,Mj,Mv,nDOF,nDOF);