function [X]=genaration_matrix_Gn_binary(U,N)
n=log2(N);
F=[1 0;1 1];
Fn=F;
for i=1:n-1
    Fn=kron(Fn,F);
end
G = Fn;
X = mod(U*G,2);
% X = bitrevorder(X);