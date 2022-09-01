function [TrigsInterpol] = SubtrigInterpol(n)
Nfactors      = (n+1)*(n+2)/2;
TrigsInterpol = zeros(Nfactors,3);
for j = 0:n
    list  = (0:j)'/n;
    index = (j*(j+1)/2+1:(j+1)*(j+2)/2)';
    TrigsInterpol(index,2) = list(end:-1:1);
    TrigsInterpol(index,3) = list;
end
TrigsInterpol(:,1 ) = 1-TrigsInterpol(:,2)-TrigsInterpol(:,3);