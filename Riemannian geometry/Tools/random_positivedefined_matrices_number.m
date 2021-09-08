function array = random_positivedefined_matrices_number(num_dimen,num_matri)

    array = zeros(num_dimen,num_dimen,num_matri);
    
    for i= 1:num_matri
        F = random_positivedefined_matrices(num_dimen);
        array(:,:,i) = F;
    end
    
end