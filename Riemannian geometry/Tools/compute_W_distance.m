function out_dis = compute_W_distance(X,Y)

        tmp_dis = trace(X + Y)-2*trace((X^(1/2)*Y*(X^(1/2)))^(1/2)); 
        out_dis = sqrt(tmp_dis);
    
end