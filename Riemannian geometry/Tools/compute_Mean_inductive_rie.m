function [out_mean,out_it] = compute_Mean_inductive_rie(X,eps)% eps is precision of algorithm

    % put matrices into array
    szx =size(X);
    num_matri = szx(3);

    % compute inital mean by geodesic of riemannian method
    mean_int = X(:,:,1)^(1/2)*((X(:,:,1)^(-1/2)*X(:,:,2)*X(:,:,1)^(-1/2))^(1/2))*(X(:,:,1)^(1/2));
    k = 2;
    i = 1;
    o = 3;
    dist_m = 10;
    
    % iteration of the mean by geodesic of riemannian method
    while dist_m > eps

        k = k + 1;
        mean_k = mean_int;
        mat = X(:,:,o);
        mean_int = mean_int^(1/2)*((mean_int^(-1/2)*mat*mean_int^(-1/2))^(1/k))*(mean_int^(1/2));
        mean_k1 = mean_int;
        dist_m = compute_riemannian_distance(mean_k,mean_k1,'A');
        if o == num_matri
            o = 1;
        else
            o = o+1;
        end
        i = i + 1;

    end
    out_it =i; %not necessary, for showing the number of iterations
    out_mean = mean_k1; % output mean
    
end