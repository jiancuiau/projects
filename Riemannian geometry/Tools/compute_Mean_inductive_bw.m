function [out_mean,out_it] = compute_Mean_inductive_bw(X,eps) % X is an array containing matrices; eps is precision of algorithm

    % put matrices into array
    szx =size(X);
    num_matri = szx(3);
    matrixtmp1 = X(:,:,1)*X(:,:,2);
    matrixtmp1 = modify_negative_matrix(matrixtmp1);
    
    matrixtmp2 = X(:,:,2)*X(:,:,1);
    matrixtmp2 = modify_negative_matrix(matrixtmp2);
    
    
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
        
        matrixtmp3 = mean_int*mat;
        matrixtmp3 = modify_negative_matrix(matrixtmp3);

        matrixtmp4 = mat*mean_int;
        matrixtmp4 = modify_negative_matrix(matrixtmp4);

        mean_int =  (1-1/k)^2*mean_int+(1/k)^2*mat+1/k*(1-1/k)*((matrixtmp3)^(1/2)+(matrixtmp4)^(1/2));
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