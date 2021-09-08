function out_dis = compute_riemannian_distance(X,Y,type)
    switch type
        case 'A'    % Affine Invariant Riemannian Metric [1,2]
            tmpEig =  eig(X,Y);  
            tmp_dis = sum(log(tmpEig).^2); 
            
        case 'S'    % Stein divergence [1,3]
            t = log(det(0.5*(X+Y))) - 0.5* (log(det(X)) + log(det(Y)));
            if t == inf || isnan(t)
                eigX = eig(X);
                eigX (eigX <= 0) = eps;
                eigY = eig(Y);
                eigY (eigY <= 0) = eps;
                t = real(sum(log(eig(0.5*(X+Y))))) - 0.5 *(real(sum(log(eigX))) + real(sum(log(eigY))));    
            end
            tmp_dis = t;
            
        case 'J'    % Jeffrey divergence [1,4]
            tmp_dis = 0.5*trace(inv(X)*Y)+0.5*trace(inv(Y)*X) - size(X,1);
            
        case 'L'    % Log-Euclidean Metric [1,5]
            tmp_dis = norm((logm(X)-logm(Y)),'fro').^2;
    end
    if tmp_dis <= 1e-15
        out_dis = 0;
    else
        out_dis = sqrt(tmp_dis);
    end
    
end