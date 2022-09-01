function [sum] = Addkron(A,B)
sum = kron(A,0*B+1)+kron(0*A+1,B);
end