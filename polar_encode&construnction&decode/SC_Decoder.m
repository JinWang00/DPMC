function U_h = SC_Decoder(N,LLR,frozenflag,UF)
n = log2(N);
L = zeros(n+1, N);        % llr
state = zeros(1, 2*N-1);  
u = zeros(n+1, N);        
U_h = zeros(1, N);
L(1,:) = LLR;             
node = 0;                 
depth = 0;                
done = 0;                 
FPos = 1;            

f = @(a,b) sign(a).*sign(b).*min(abs(a),abs(b)); % min-sum  
g = @(a,b,c) (1-2*c).*a + b; 

while(done == 0)
    if depth == n 
        if frozenflag(node+1) == 1     % forzen bits
            u(n+1, node+1) = UF(FPos); 
            FPos = FPos + 1;
        else
            if L(n+1, node+1) >= 0     
                u(n+1, node+1) = 0;
            else
                u(n+1, node+1) = 1;
            end            
        end
        if node == (N-1) 
            U_h = u(n+1,:);
            %U_h = bitrevorder(U_h); % bit reverse
            done = 1;
        else
            node = floor(node/2); 
            depth = depth - 1;    
        end
    else 
        npos = 2^depth + node;  
        if state(npos) == 0     
            temp = 2^(n-depth); 
            Ln = L(depth+1, temp*node+1:temp*(node+1));
            a = Ln(1:temp/2);
            b = Ln(temp/2+1:end); 
            temp = temp/2;
            node = node*2; 
            depth = depth + 1;
            L(depth+1, temp*node+1:temp*(node+1)) = f(a,b); 
            state(npos) = 1;
        else
            if state(npos) == 1  
                temp = 2^(n-depth); 
                Ln = L(depth+1, temp*node+1:temp*(node+1));
                a = Ln(1:temp/2);
                b = Ln(temp/2+1:end);
                lnode = 2*node; 
                ldepth = depth + 1;
                ltemp = temp/2;
                un = u(ldepth+1, ltemp*lnode+1:ltemp*(lnode+1)); 
                temp = temp/2;
                node = node*2+1;
                depth = depth + 1;
                L(depth+1, temp*node+1:temp*(node+1)) = g(a,b,un); 
                state(npos) = 2;
            else 
                temp = 2^(n-depth);
                lnode = 2*node;    
                rnode = 2*node + 1; 
                ndepth = depth + 1;
                ntemp = temp/2;    
                ul = u(ndepth+1, ntemp*lnode+1:ntemp*(lnode+1)); 
                ur = u(ndepth+1, ntemp*rnode+1:ntemp*(rnode+1));
                u(depth+1, temp*node+1:temp*(node+1)) = [mod(ul + ur, 2) ur];
                node = floor(node/2);
                depth = depth - 1;
            end
        end
    end
end
