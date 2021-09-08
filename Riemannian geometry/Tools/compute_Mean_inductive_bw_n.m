function [out_mean,out_it] = compute_Mean_inductive_bw_n(X,eps) % X is an array containing matrices; eps is precision of algorithm

    % put matrices into array
    szx =size(X);
    num_matri = szx(3);
    matrixtmp1 = X(:,:,1)*X(:,:,2);
    matrixtmp1 = modify_negative_matrix(matrixtmp1);
    
    matrixtmp2 = X(:,:,2)*X(:,:,1);
    matrixtmp2 = modify_negative_matrix(matrixtmp2);
    inv_mat = zeros(szx);
    for i = 1:num_matri
    inv_mat(:,:,i)= X(:,:,i)^(-1);
    end
    % compute inital mean by geodesic of bw method
    mean_int =  (1-1/2)^2*X(:,:,1)+(1/2)^2*X(:,:,2)+1/2*(1-1/2)*((matrixtmp1)^(1/2)+(matrixtmp2)^(1/2));

    k = 2;
    i = 1;
    o = 3;
    dist_m = 10;
    % the iteration of means by geodesic of bw method
    while dist_m > eps

        k = k + 1;
        mean_k = mean_int;
        mat = X(:,:,o);
        
        matrixtmp3 = mat*mean_int;
        matrixtmp3 = modify_negative_matrix(matrixtmp3);
        sq_matrixtmp3 = (matrixtmp3)^(1/2);
        
        sq_matrixtmp4 = inv_mat(:,:,o)*sq_matrixtmp3*mat;

        mean_int =  (1-1/k)^2*mean_int+(1/k)^2*mat+1/k*(1-1/k)*(sq_matrixtmp4+sq_matrixtmp3);
        mean_k1 = mean_int;
        dist_m = compute_W_distance(mean_k,mean_k1);
        
        % computing the number of iteration
        if o == num_matri 
            o = 1;
        else
            o = o+1;
        end
        i = i + 1;

    end
    out_it =i; % not necessary, computing the number of iteration
    out_mean = mean_k1; % mean for output
    
end