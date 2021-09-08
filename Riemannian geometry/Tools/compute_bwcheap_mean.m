function center = compute_bwcheap_mean(matrices)

szx =size(matrices);
num_matri = szx(3);
num_dimens = szx(1);

meanarray = zeros(num_dimens,num_dimens,num_matri);

iter = 0;

dista = ones(num_matri,1);

    while all(dista(:) > 1e-6)% 1e-6 is the maximum elements' difference among means

        iter = iter +1 ;
        % compute means by every matrix's tangent space
        for l = 1:num_matri

            M = matrices(:,:,l);
            [~,num_de,num_mat] = size(matrices);
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

            for i = 1:1:num_de
                for r = 1:1:num_de
                    va =(eigv(i)+eigv(r))^(-1);
                    W(i,r) = va;
                end
            end
            
            % projection to manifold
            M = V*(W.*(V_1*S*V+2*D))*D*(W.*(V_1*S*V+2*D))*V_1;

            meanarray(:,:,l) = M;

        end
        
         % determine convergence by the maximum element's difference
        
        center =  mean(meanarray,3);
        
        for i = 1:num_mat
            matabs = abs(matrices(:,:,1)-center);
            matmax = max(matabs(:));
            dista(i,1)= matmax;
        end

        matrices =  meanarray;

    end


end