function out_sigma1 = random_positivedefined_matrices(num_dimen)
    A=randn(num_dimen);
    A=A'*A;
    B=0.01*eye(num_dimen);
    out_sigma1 = A+B;
end
    