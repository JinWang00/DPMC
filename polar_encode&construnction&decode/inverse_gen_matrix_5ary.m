function [X]=inverse_gen_matrix_5ary(U,N)
n=log2(N);
F=[1 0;4 1]; % F=[1 0;-1 1]
Fn=F;
for i=1:n-1
    Fn=kron(Fn,F);
end
G = Fn;
X = mod(U*G,5);
% X = bitrevorder(X);