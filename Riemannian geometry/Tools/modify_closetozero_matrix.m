function matrix = modify_closetozero_matrix(mat)
    sez = size(mat);
    eig1 = ones(sez(1),1);
    eig1 = -eig1;
    vec1 = ones(sez(1),1);
    vec1 = 1e-14*vec1;
    mat1 =diag(vec1);
    while all(eig1 > 1e-14) == false 
                [vec1,eig1]= eig(mat);
                 eig1 = abs(eig1+mat1);
                 mat = vec1*eig1/vec1;
                  matrix = mat;
                 eig1 = eig(matrix);
    end
end