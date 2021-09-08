function B = eigvalues_closeto_zero_matrix(num_dimen)


    A = random_positivedefined_matrices(num_dimen);


    [V,D] = eig(A);

    D = diag(D);
        
    eps = 1e-7; % the level of eigen values 
        
    D(num_dimen) = D(num_dimen)*eps;

    D = diag(D);

    B = V*D/V;

end