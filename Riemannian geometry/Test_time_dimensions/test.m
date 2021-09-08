
t = 1/4;

mean1 = (1-t)^2*T2+t^2*T3+t*(1-t)*((T2*T3)^(1/2)+(T3*T2)^(1/2));

dis1 = compute_W_distance(mean1,T2); % T2 closer
dis2 = compute_W_distance(mean1,T3);


%%
num_dimens = 10;
num_matri = 10;
esp = 1e-2;

data1 = random_positivedefined_matrices_number(num_dimens,num_matri);

tic
mean_3 = compute_Mean_inductive_bw(data1,esp);
tim_m1 = toc;
    
tic
mean_4 = compute_Mean_inductive_rie(data1,esp);
tim_m2 = toc;


%% function cheap bw


num_dimens = 10;
num_matri = 10;

matrices = random_positivedefined_matrices_number(num_dimens,num_matri);

meanarray = zeros(num_dimens,num_dimens,num_matri);

iter = 0;

dista = ones(num_matri,1);

while all(dista(:) > 1e-6) 
    iter = iter +1 ;
    
    for l = 1:num_matri

        M = matrices(:,:,l);
        eps = 1e-6;
        [~,num_de,num_mat] = size(matrices);

        M_o = M;
        [V,D] = eig(M);
        V_1 = V^(-1);
        eigv = eig(M);
        s_matrix = zeros(num_de,num_de,num_mat);

        for i = 1:1:num_mat
            p = matrices(:,:,i);
            s_matrix(:,:,i) = (M * p)^(1/2)+(p * M)^(1/2) - 2 * M;
        end

        S = mean(s_matrix,3);

        W = zeros(num_de,num_de);

        for i = 1:1:num_de
            for r = 1:1:num_de
                va =(eigv(i)+eigv(r))^(-1);
                W(i,r) = va;
            end
        end

        M = V*(W.*(V_1*S*V+2*D))*D*(W.*(V_1*S*V+2*D))*V_1;

        meanarray(:,:,l) = M;


    end
        center =  mean(meanarray,3);
        for i = 1:num_mat
            matabs = abs(matrices(:,:,1)-center);
            matmax = max(matabs(:));
            dista(i,1)= matmax;
        end

    matrices =  meanarray;

end 



%% cheap rie test 
num_dimens = 10;
num_matri = 10;

matrices = random_positivedefined_matrices_number(num_dimens,num_matri);



meanarray = zeros(num_dimens,num_dimens,num_matri);

iter = 0;

dista = ones(num_matri,1);

distable = ones(num_matri,100);

while all(dista(:) > 1e-6)
    iter = iter +1 ;
    
   for l = 1:num_matri

         M = matrices(:,:,l);
        A = M ^ (1/2);  
        B = A ^ (-1);       

        S = zeros(size(M));
        for j_th = 1 : num_matri
            C = matrices(:,:,j_th);
            S = S + A * logm(B * C * B) * A;
        end
        S = S / num_matri;

        M = A * expm(B * S * B) * A; 

        meanarray(:,:,l) = M;

    end

    center =  mean(meanarray,3);
    for i = 1:num_mat
        matabs = abs(matrices(:,:,1)-center);
        matmax = max(matabs(:));
        dista(i,1)= matmax;
    end
    
    distable(:,iter)= dista;
    matrices =  meanarray;
    
end    
    
%% test bw and rie function
% number of dimensions of matrices
num_dimens = 20;
% number of matrices
num_matri = 10;

% generate symmetric positive definite matrices
matrices = random_positivedefined_matrices_number(num_dimens,num_matri);

% compute means by different methods 

mean1 = compute_bwcheap_mean(matrices);
mean2 = compute_riecheap_mean(matrices);

mean3 = compute_Mean_inductive_bw(matrices,1e-2);

tic
mean4 = compute_Mean_inductive_rie(matrices,1e-2);
toc
toc3=toc;
mean5 = compute_bwpro_mean(matrices);
mean6 = compute_riepro_mean(matrices,'A',1000);

mean7 = compute_Mean_inductive_bw_n(matrices,1e-2);

tic
mean8 = compute_Mean_inductive_rie_n(matrices,1e-2);
toc
toc4 = toc;

%% 
Sn = matrices(:,:,1);
A = matrices(:,:,2);
D = Sn^(1/2)*((Sn^(-1/2)*A*Sn^(-1/2))^(1/3))*(Sn^(1/2));
Sd = compute_riemannian_distance(Sn,D,'A');
Ad = compute_riemannian_distance(A,D,'A');

%%
Sn = matrices(:,:,1);
A = matrices(:,:,2);

H1 = (Sn*A)^(1/2);
H2 = A^(-1)*(A*Sn)^(1/2)*A;
%%
M = data1(:,:,1);
p = data1(:,:,2);
test = (M*p)^(1/2)
test1 = sqrtm(M*p)