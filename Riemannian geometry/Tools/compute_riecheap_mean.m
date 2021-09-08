function center = compute_riecheap_mean(matrices)

szx =size(matrices);
num_matri = szx(3);
num_dimens = szx(1);

meanarray = zeros(num_dimens,num_dimens,num_matri);

iter = 0;

dista = ones(num_matri,1);

distable = ones(num_matri,100);

    while all(dista(:) > 1e-6)% 1e-6 is the maximum element-wised difference among means
        iter = iter +1 ;
       % compute means by every matrix's tangent space
           for l = 1:num_matri

                 M = matrices(:,:,l);
                A = M ^ (1/2);  
                B = A ^ (-1);       
                S = zeros(size(M));
                % projection to tangent space
                for j_th = 1 : num_matri
                    C = matrices(:,:,j_th);
                    S = S + A * logm(B * C * B) * A;
                end
                S = S / num_matri;
                % projection to manifold
                M = A * expm(B * S * B) * A; 

                meanarray(:,:,l) = M;

           end
        % determine convergence by the maximum elements' difference
        center =  mean(meanarray,3);

        for i = 1:num_matri
            matabs = abs(matrices(:,:,1)-center);
            matmax = max(matabs(:));
            dista(i,1)= matmax;
        end

        distable(:,iter)= dista; % not not necessary, for checking process.

        matrices =  meanarray;

    end    
    
end