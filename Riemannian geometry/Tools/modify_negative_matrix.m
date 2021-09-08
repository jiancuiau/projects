function matrix = modify_negative_matrix(mat)
    sez = size(mat);
    eig1 = ones(sez(1),1);
    eig1 = -eig1;
    while all(eig1 > 0) == false 
                [vec1,eig1]= eig(mat);
                 eig1 = abs(eig1);
                 mat = vec1*eig1/vec1;
                  matrix = mat;
                 eig1 = eig(matrix);
    end
end