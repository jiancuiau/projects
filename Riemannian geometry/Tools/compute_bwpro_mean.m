function M = compute_bwpro_mean(matrices)


% compute inital mean
    M = mean(matrices,3);
    
% precision of algorithm
    eps = 1e-6;
    
    [~,num_de,num_mat] = size(matrices);
    dis = 10;
    n = 0;
    
    while dis > eps

        n = n + 1;
        M_o = M;
        [V,D] = eig(M);
        V_1 = V^(-1);
        eigv = eig(M);
        s_matrix = zeros(num_de,num_de,num_mat);
        
        % projection to tangent space
        for i = 1:1:num_mat
            p = matrices(:,:,i);
         
            s_matrix(:,:,i) = (M * p)^(1/2)+(p * M)^(1/2) - 2 * M;
        end
        
        S = mean(s_matrix,3);

        W = zeros(num_de,num_de);

        for l = 1:1:num_de
            for r = 1:1:num_de
                va =(eigv(l)+eigv(r))^(-1);
                W(l,r) = va;
            end
        end
        
        % projection to manifold 
        M = V*(W.*(V_1*S*V+2*D))*D*(W.*(V_1*S*V+2*D))*V_1;

        dis = compute_W_distance(M,M_o);
    
     end
    
end