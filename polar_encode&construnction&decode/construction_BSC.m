function [PeLast,Indices] = construction_BSC(N,p)

PeLast = zeros(1,N);
colmax = 32; 
n = log2(N);
Pmat = struct('mat',cell(1,N)); 
Pmat(1).mat = [1;p];

for K = 1:n
    for I = 1:2^(K-1)
        temp1 = Pmat(I).mat(1,:).*Pmat(I).mat(2,:);
        temp = [Pmat(I).mat(1,:)-temp1 ; temp1];
        temp = kron(temp,temp);
        Pup = [temp(1,:)+temp(4,:);temp(3,:)+temp(2,:)];
        Pdown = [temp(1,:),temp(3,:);temp(4,:),temp(2,:)];
        Pu = [];
        Pd = [];
        Pu(1,:) = sum(Pup);
        Pu(2,:) = min(Pup)./Pu(1,:);
        Pd(1,:) = sum(Pdown);
        Pd(2,:) = min(Pdown)./Pd(1,:);
        
        % Pu
        % merge same col
        idx = Pu(1,:)==0;
        Pu(:,idx) = [];
        temp = unique(Pu(2,:));
        col = length(temp);
        Tu = zeros(2,col);
        for J = 1:col
            idx = (Pu(2,:) == temp(J));
            Tu(1,J) = sum(Pu(1,idx));
            Tu(2,J) = temp(J);
        end

        while col > colmax
            minDiffC = 1; 
            minIdx = 0;
            for J = 1:col-1 
                t1 = Tu(1,J);
                t2 = Tu(1,J+1);
                t3 = Tu(2,J);
                t4 = Tu(2,J+1);
                C1 = (t1 + t2)*getEntropy( (t1*t3 + t1*t4) / (t1+t2) );
                C2 = t1*getEntropy(t3) + t2*getEntropy(t4);
                
                if C1-C2 < minDiffC %find minimum delta I
                    minDiffC = C1-C2; 
                    minIdx = J;
                end
            end
            Tu(2,minIdx) = ((Tu(1,minIdx)*Tu(2,minIdx)+Tu(1,minIdx+1)*Tu(2,minIdx+1)) / (Tu(1,minIdx)+Tu(1,minIdx+1)));
            Tu(1,minIdx) = (Tu(1,minIdx)+Tu(1,minIdx+1));
            Tu(:,minIdx+1) = [];
            col = col-1;
        end
        Pmat(I).mat = Tu;
        
        idx = Pd(1,:)==0;
        Pd(:,idx) = [];
        temp = unique(Pd(2,:));
        col = length(temp);
        Td = zeros(2,col);
        for J = 1:col
            idx = (Pd(2,:) == temp(J));
            Td(1,J) = sum(Pd(1,idx));
            Td(2,J) = temp(J);
        end
        
        %merge
        while col > colmax
            minDiffC = 1;
            minIdx = 0;
            for J = 1:col-1 
                C1 = (Td(1,J)+Td(1,J+1))*getEntropy( (Td(1,J)*Td(2,J)+Td(1,J+1)*Td(2,J+1)) / (Td(1,J)+Td(1,J+1)) );
                C2 = Td(1,J)*getEntropy(Td(2,J)) + Td(1,J+1)*getEntropy(Td(2,J+1));
                if C1-C2 < minDiffC
                    minDiffC = C1-C2;
                    minIdx = J;
                end
            end
            Td(2,minIdx) = ((Td(1,minIdx)*Td(2,minIdx)+Td(1,minIdx+1)*Td(2,minIdx+1)) / (Td(1,minIdx)+Td(1,minIdx+1)));
            Td(1,minIdx) = (Td(1,minIdx)+Td(1,minIdx+1));
            Td(:,minIdx+1) = [];
            col = col-1;
        end
        Pmat(I+2^(K-1)).mat = Td;
    end
end

for I = 1:N
    PeLast(I) = sum((Pmat(I).mat(1,:)).*Pmat(I).mat(2,:));
end
idx = bin2dec(fliplr(dec2bin((1:N)-1)))+1;
PeLast = PeLast(idx);
[~,Indices] = sort(PeLast,2);












