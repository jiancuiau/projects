function out = create_SPD(sigma1,num_dimens,num_matri)

eig(sigma1);
df = 2*num_dimens;
[~,D] = wishrnd(sigma1,df);
data1 = zeros(num_dimens,num_dimens,num_matri);


    for i= 1:num_matri
        F = wishrnd(sigma1,df,D)/df;
        data1(:,:,i) = F;
    end
    
out = data1;
    