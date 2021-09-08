function array = noise_matrix_SPD(num_dimens,num_matrix)

array = zeros(num_dimens,num_dimens,num_matrix);


    for i = 1:num_matrix

       tmp = random_positivedefined_matrices(num_dimens);
        
        array(:,:,i) = tmp;

    end

end
