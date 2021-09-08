function array = noise_matrix(num_dimens,num_matrix)

array = zeros(num_dimens,num_dimens,num_matrix);


    for i = 1:num_matrix

        mu=0;
        sigma=1;
        
        index = reshape(triu(ones(num_dimens)),num_dimens*num_dimens,1)==1;

        vecback = zeros(1,length(index));

        tt = num_dimens*(num_dimens+1)/2;

        n = 0;

        out = sigma *randn(1,tt)+mu;


        for a = 1:length(index)
            if index(a) == 1
                n = n + 1;
                vecback(a) = out(n);
            else
                 vecback(a) = 0;
            end
        end
        tmp = reshape(vecback,[],num_dimens);
        tmp2 = triu(tmp,1) +diag(sqrt((diag(tmp)).^2))+ tril(tmp.',-1);
        
        array(:,:,i) = tmp2;

    end

end
