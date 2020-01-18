%% Linear Discrite Time System Matrices
function [H,L,M,Z] = LDTS_MPC_matrices(A,b,c,Ku,Ky,lamda)
%% L matrices
L = 2*eye(Ku+1);
L(end,end) = 1;
for i=1:Ku
    for j=1:Ku
        L(i+1,i) = -1;
        L(i,i+1) = -1;
    end
end
if Ku == 0
    L = 2;
end

%% Z matrices
Z = [];
for k=1:Ky
    Zn = c'*A^k;
    Z = [Z;Zn];
end
%% M matrices
M = [c'*b,zeros(1,Ku)];
for i=2:Ky
    M(i,:) = [c'*A^(i-1)*b,M(i-1,1:end-1)];
end
for k = Ku+2:Ky
    % M(k,end) = M(k-1,end-1) + M(k-1,end);
    aa = 0;
    for l=1:k-Ku-1
        aa = aa + c'*A^(l-1)*b;
    end
    M(k,Ku+1) = M(k,Ku+1) + aa;
end
%% H matrices
H = 2*(M'*M + lamda * L );



end

