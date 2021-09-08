%%

clc;

clear;

num_dimens = 3; % number of dimensions

matrix = eigvalues_closeto_zero_matrix(num_dimens);

num_mat = 100; % number of matrices

mat_array = zeros(num_dimens,num_dimens,num_mat);

noise = [0.00001,0,0; % noise matrice 
        0,0,0;
        0,0,0];

% create positive-defined matrices

for i= 1:100
    
    isposdef = 0;

    while isposdef == 0

        matrix = matrix + noise;

        d = eig(matrix);

        isposdef = all(d > 0);

    end

    mat_array(:,:,i) = matrix;

end

T1 = mat_array(:,:,1);

T2 = random_positivedefined_matrices(num_dimens);

T3 = random_positivedefined_matrices(num_dimens);

matrices = zeros(num_dimens,num_dimens,3);

matrices(:,:,1) =  T1;
matrices(:,:,2) =  T2;
matrices(:,:,3) =  T3;



%% projection algorithm

r_m = compute_riepro_mean(matrices,'A',10000);

bw_m = compute_bwpro_mean(matrices);

table = zeros(num_mat-1,3);

meanarray_r = zeros(num_dimens,num_dimens,num_mat-1);

meanarray_r(:,:,1) = r_m;

meanarray_bw = zeros(num_dimens,num_dimens,num_mat-1);

meanarray_bw(:,:,1) = bw_m;

% compute distance between means

for i = 2:num_mat

    T1 = mat_array(:,:,i);

    matrices(:,:,1) =  T1;


    r_m_t = compute_riepro_mean(matrices,'A',10000);
    
    meanarray_r(:,:,i) = r_m_t;

    bw_m_t = compute_bwpro_mean(matrices);
    
    meanarray_bw(:,:,i) = bw_m_t;

    table(i-1,1)= i-1;

    table(i-1,2)= compute_W_distance(r_m_t,r_m);

    table(i-1,3)= compute_W_distance(bw_m_t,bw_m);

    r_m = r_m_t;

    bw_m = bw_m_t;

end


x = table(:,1);

a = table(:,2);

b = table(:,3);


figure(1);  
plot(x,a,'-*b', 'linewidth', 1.1)
set(gca,'FontSize',20);
xlabel('Number','fontsize',20)  
ylabel('distance of rie mean','fontsize',20)  

figure(2);  
plot(x,b,'-*b', 'linewidth', 1.1)
set(gca,'FontSize',20);
xlabel('Number','fontsize',20)  
ylabel('distance of bw mean','fontsize',20) 


%% cheap2 algorithm

r_m = compute_riecheap_mean(matrices);

bw_m = compute_bwcheap_mean(matrices);

table = zeros(num_mat-1,3);

meanarray_r = zeros(num_dimens,num_dimens,num_mat-1);

meanarray_r(:,:,1) = r_m;

meanarray_bw = zeros(num_dimens,num_dimens,num_mat-1);

meanarray_bw(:,:,1) = bw_m;

% compute distance between means

for i = 2:num_mat
    
    i

    T1 = mat_array(:,:,i);

    matrices(:,:,1) =  T1;


    r_m_t = compute_riecheap_mean(matrices);
    
    meanarray_r(:,:,i) = r_m_t;

    bw_m_t = compute_bwcheap_mean(matrices);
    
    meanarray_bw(:,:,i) = bw_m_t;

    table(i-1,1)= i-1;

    table(i-1,2)= rie_dis(r_m_t,r_m);

    table(i-1,3)= compute_W_distance(bw_m_t,bw_m);

    r_m = r_m_t;

    bw_m = bw_m_t;

end


x = table(:,1);

a = table(:,2);

b = table(:,3);


figure(1);  
plot(x,a,'-*b', 'linewidth', 1.1)
set(gca,'FontSize',20);
xlabel('Number','fontsize',20)  
ylabel('distance of rie mean','fontsize',20)  

figure(2);  
plot(x,b,'-*b', 'linewidth', 1.1)
set(gca,'FontSize',20);
xlabel('Number','fontsize',20)  
ylabel('distance of bw mean','fontsize',20) 

%% inductive algorithm
eps = 1e-2;

r_m = compute_Mean_inductive_rie(matrices,eps);

bw_m = compute_Mean_inductive_bw(matrices,eps);

table = zeros(num_mat-1,3);

meanarray_r = zeros(num_dimens,num_dimens,num_mat-1);

meanarray_r(:,:,1) = r_m;

meanarray_bw = zeros(num_dimens,num_dimens,num_mat-1);

meanarray_bw(:,:,1) = bw_m;

% compute distance between means

for i = 2:num_mat

    T1 = mat_array(:,:,i);

    matrices(:,:,1) =  T1;


    r_m_t = compute_Mean_inductive_rie(matrices,eps);
    
    meanarray_r(:,:,i) = r_m_t;

    bw_m_t = compute_Mean_inductive_bw(matrices,eps);
    
    meanarray_bw(:,:,i) = bw_m_t;

    table(i-1,1)= i-1;

    table(i-1,2)= compute_W_distance(r_m_t,r_m);

    table(i-1,3)= compute_W_distance(bw_m_t,bw_m);

    r_m = r_m_t;

    bw_m = bw_m_t;

end


x = table(:,1);

a = table(:,2);

b = table(:,3);


figure(1);  
plot(x,a,'-*b', 'linewidth', 1.1)
set(gca,'FontSize',20);
xlabel('Number','fontsize',20)  
ylabel('distance of rie mean','fontsize',20)  

figure(2);  
plot(x,b,'-*b', 'linewidth', 1.1)
set(gca,'FontSize',20);
xlabel('Number','fontsize',20)  
ylabel('distance of bw mean','fontsize',20) 





%% cheap algorithm

r_m = compute_riecheap_mean(matrices);

bw_m = compute_bwcheap_mean(matrices);

table = zeros(num_mat-1,3);

meanarray_r = zeros(num_dimens,num_dimens,num_mat-1);

meanarray_r(:,:,1) = r_m;

meanarray_bw = zeros(num_dimens,num_dimens,num_mat-1);

meanarray_bw(:,:,1) = bw_m;

% compute distance between means

for i = 2:num_mat

    T1 = mat_array(:,:,i);

    matrices(:,:,1) =  T1;


    r_m_t = compute_riecheap_mean(matrices);
    
    meanarray_r(:,:,i) = r_m_t;

    bw_m_t = compute_bwcheap_mean(matrices);
    
    meanarray_bw(:,:,i) = bw_m_t;

    table(i-1,1)= i-1;

    table(i-1,2)= rie_dis(r_m_t,r_m);

    table(i-1,3)= compute_W_distance(bw_m_t,bw_m);

    r_m = r_m_t;

    bw_m = bw_m_t;

end


x = table(:,1);

a = table(:,2);

b = table(:,3);


figure(1);  
plot(x,a,'-*b', 'linewidth', 1.1)
set(gca,'FontSize',20);
xlabel('Number','fontsize',20)  
ylabel('distance of rie mean','fontsize',20)  

figure(2);  
plot(x,b,'-*b', 'linewidth', 1.1)
set(gca,'FontSize',20);
xlabel('Number','fontsize',20)  
ylabel('distance of bw mean','fontsize',20) 

%%
test_mat(:,:,1) = T1; 
test_mat(:,:,2) = T2; 
test_mat(:,:,3) = T3; 

compute_cheap_mean(test_mat);

matrices = test_mat;

M = matrices(:,:,1);
    eps = 1e-6;
    [~,num_de,num_mat] = size(matrices);
    dis = 10;
    n = 0;
    
    while dis > eps

        n = n + 1;
        M_o = M;
        [V,D] = eig(M);
        
        V_1 = V^(-1);
        eigv = eig(M);
        s_matrix = zeros(num_de,num_de,num_mat);
        
        for i = 1:1:num_mat
            p = matrices(:,:,i);
            matrixtmp1 = M * p;
             matrixtmp1 = modify_negative_matrix(matrixtmp1);
            matrixtmp2 = p * M;
             matrixtmp2 = modify_negative_matrix(matrixtmp2);
            s_matrix(:,:,i) = (matrixtmp1)^(1/2)+(matrixtmp2)^(1/2) - 2 * M;
        end
        %if isreal(s_matrix) == false
               % break
        %end
        S = mean(s_matrix,3);

        W = zeros(num_de,num_de);

        for l = 1:1:num_de
            for r = 1:1:num_de
                va =(eigv(l)+eigv(r))^(-1);
                W(l,r) = va;
            end
        end

        M = V*(W.*(V_1*S*V+2*D))*D*(W.*(V_1*S*V+2*D))*V_1;
        %if isreal(M) == false
                %break
        %end
        
         
        dis = compute_W_distance(M,M_o);

    end
%% test rie
tablerieeig = zeros(99,3);
tablerietest = zeros(99,2);
for i =1:99
    X = meanarray_r(:,:,i);
    Y = meanarray_r(:,:,i+1);
    tmpEig =  eig(X,Y); 
    tablerieeig(i,1:3) = tmpEig;
    if isreal(tmpEig) == false
                break
    end
    %if all(tmpEig > 0) == false
     %   break
    %end
    tmp_dis = sum(log(tmpEig).^2);
    tablerietest(i,1) = tmp_dis;
    out =sqrt(tmp_dis);
    tablerietest(i,2) = out;
    %if isreal(out) == false
     %           break
    %end
end

%% 

 eig1 = eig(T1);
 eig2 = eig(T1^2);
 eig3 = eig(T1*T2);
 eig9 = eig(T2*T1);
 eig10 = eig(T1*T3);
 eig11 = eig(T3*T1);
 
 eig4 = eig(T2);
 eig5 = eig(T3);
 eig6 = eig(T2*T3);
 eig7 = eig((T2*T3)^(1/2));
 eig8 = eig((T2)^(1/2)*(T3)^(1/2));
 
 
%%

eig(T2);
eig(M);
eig(T2*M);%% error
eig(M*T2);

test_matrix1= T2;
test_matrix2= M;
eig(test_matrix1)
eig(test_matrix2)
eig(test_matrix1*test_matrix2)%% error with 1e-20 matrix
eig(test_matrix2*test_matrix1)



%%
r_m = rie_mean(matrices,10000);

bw_m = compute_cheap_mean(matrices);

table = zeros(num_mat-1,3);


meanarray_bw = zeros(num_dimens,num_dimens,num_mat-1);

meanarray_bw(:,:,1) = bw_m;


 T1 = mat_array(:,:,2);

 matrices(:,:,1) =  T1;


    bw_m_t = compute_bwpro_mean(matrices);
    
    meanarray_bw(:,:,2) = bw_m_t;
    
    
%% test for rie real results

%%
f = 0;
tims = 10;
for j = 1:tims

num_dimens = 3; % number of dimensions

matrix = eigvalues_closeto_zero_matrix(num_dimens);

eigs = eig(matrix);

ratioeig = max(eigs)/min(eigs);

num_mat = 5; % number of matrices

mat_array = zeros(num_dimens,num_dimens,num_mat);

noise = [0.00001,0,0; % noise matrice 
        0,0,0;
        0,0,0];

% create positive-defined matrices

for i= 1:num_mat
    
    isposdef = 0;

    while isposdef == 0

        matrix = matrix + noise;

        d = eig(matrix);

        isposdef = all(d > 0);

    end

    mat_array(:,:,i) = matrix;

end

T1 = mat_array(:,:,1);

T2 = random_positivedefined_matrices(num_dimens);

T3 = random_positivedefined_matrices(num_dimens);

matrices = zeros(num_dimens,num_dimens,3);

matrices(:,:,1) =  T1;
matrices(:,:,2) =  T2;
matrices(:,:,3) =  T3;

i


r_m = compute_riepro_mean(matrices,'A',10000);

bw_m = compute_bwpro_mean(matrices);

table = zeros(num_mat-1,3);

meanarray_r = zeros(num_dimens,num_dimens,num_mat-1);

meanarray_r(:,:,1) = r_m;

meanarray_bw = zeros(num_dimens,num_dimens,num_mat-1);

meanarray_bw(:,:,1) = bw_m;

% compute distance between means

for i = 2:num_mat

    T1 = mat_array(:,:,i);

    matrices(:,:,1) =  T1;


    r_m_t = compute_riepro_mean(matrices,'A',10000);
    
    meanarray_r(:,:,i) = r_m_t;

    bw_m_t = compute_bwpro_mean(matrices);
    
    meanarray_bw(:,:,i) = bw_m_t;

    table(i-1,1)= i-1;

    table(i-1,2)= rie_dis(r_m_t,r_m);

    table(i-1,3)= compute_W_distance(bw_m_t,bw_m);

    r_m = r_m_t;

    bw_m = bw_m_t;

end


x = table(:,1);

a = table(:,2);

b = table(:,3);

    if isreal(meanarray_r) == false
        f = f + 1;
    end
end 

complratio = f/tims

%%
figure(1);  
plot(x,a,'-*b', 'linewidth', 1.1)
set(gca,'FontSize',20);
xlabel('Number','fontsize',20)  
ylabel('distance of rie mean','fontsize',20)  

figure(2);  
plot(x,b,'-*b', 'linewidth', 1.1)
set(gca,'FontSize',20);
xlabel('Number','fontsize',20)  
ylabel('distance of bw mean','fontsize',20) 

%%

spd_matrices = matrices;
[dims,~,num_spd] = size(spd_matrices);
    M  = mean(spd_matrices, 3);
        max_iter = 100000;    
            for ite_th = 1 : max_iter
                A = M ^ (1/2);      %-- A = C^(1/2)
                B = A ^ (-1);       %-- B = C^(-1/2)

                S = zeros(size(M));
                for j_th = 1 : num_spd
                    C = spd_matrices(:,:,j_th);
                    D = B*C*B;
                    D = modify_negative_matrix(D);
                    E = A * logm(D) * A;
                    if isreal(E) == false
                            break
                    end
                    S = S + E;
                end
                S = S / num_spd;

                M = A * expm(B * S * B) * A; 
                    if isreal(M) == false
                            break
                    end
                eps = norm(S, 'fro');
                if (eps < 1e-6)
                    break;
                end
            end
            
%%
D = eigvalues_closeto_zero_matrix(3);
[V,U] = eig(D);
eigv = eig(D);
noisee = [1e-14,0,0;0,1e-14,0;0,0,1e-14];
eigd = U+noisee;
%V*U*V^(-1)
H = V*eigd*V^(-1);
logm(H)
 H = modify_closetozero_matrix(D);
 logm(H)
 
 %%
 
 X= meanarray_r(:,:,7);
 Y= meanarray_r(:,:,8);
 
 X = (X+X')/2;
 Y = (Y+Y')/2;
 X = modify_negative_matrix(X);
 Y = modify_negative_matrix(Y);
 tmpEig =  eig(X,Y);  
 tmp_dis = sum(log(tmpEig).^2);
            
    if tmp_dis <= 1e-15
        out_dis = 0;
    else
        out_dis = sqrt(tmp_dis);
    end
    
%%
tt = compute_riepro_mean(matrices,'A',10000);

eigt = eig(matrices(:,:,1));
%%
num_dimens = 3;
% number of matrices
num_matri = 3;

% generate symmetric positive definite matrices
matrices = random_positivedefined_matrices_number(num_dimens,num_matri);
%%
A = random_positivedefined_matrices(3);
eig(A)

%%

num_dimens = 3;

T1 = eigvalues_closeto_zero_matrix(num_dimens);

eigs = eig(T1);

ratioeig = max(eigs)/min(eigs);

T2 = random_positivedefined_matrices(num_dimens);

T3 = random_positivedefined_matrices(num_dimens);

matrices = zeros(num_dimens,num_dimens,3);

matrices(:,:,1) =  T1;
matrices(:,:,2) =  T2;
matrices(:,:,3) =  T3;



spd_matrices = matrices;
[dims,~,num_spd] = size(spd_matrices);
    M  = mean(spd_matrices, 3);
        max_iter = 100000;    
        meantable = zeros(dims,dims,1000);
        epstable = zeros(1,1000);
            for ite_th = 1 : max_iter
                A = M ^ (1/2);      %-- A = C^(1/2)
                B = A ^ (-1);       %-- B = C^(-1/2)

                S = zeros(size(M));
                for j_th = 1 : num_spd
                    C = spd_matrices(:,:,j_th);
                    D = B*C*B;
                    H = logm(D);
                    if isreal(H) == false
                        
                        F = modify_negative_matrix(D); 
                        H = logm(F);
                    end

                    if isreal(H) == false
                            break
                    end
                    S = S + A * H * A;
                end
            
                S = S / num_spd;

                M = A * expm(B * S * B) * A;
                M = (M + M.')/2;
                
                
                if isreal(M) == false
                     break
                end
                meantable(:,:,ite_th) = M;
                eps = norm(S, 'fro');
                epstable(:,ite_th) = eps;
                if (eps < 1e-6)
                    break;
                end
            end
 %%
  bw_mean = compute_bwpro_mean(matrices);
 %%
 
 [EE,EV]=  eig(D);
EVV = EV + [0,0,0; % noise matrice 
        0,0,0;
        0,0,1e-12];
EEE= EE*EVV/EE;
logm(EEE)
%%
t= 1
mattt = eigvalues_closeto_zero_matrix(3)
eigmattt = eig(mattt)
mattt = modify_negative_matrix(mattt)
logm(mattt)
eig(mattt)

%%
mattt = eigvalues_closeto_zero_matrix(3)
eigmattt = eig(mattt)
matrix = modify_negative_matrix(mattt);
eig(matrix)
logm(matrix)


%%
%% test for rie real results

%%
f = 0;
tims = 100;
rattable = zeros(tims,3);
for j = 1:tims

rattable(j,1) = j;
num_dimens = 3; % number of dimensions

matrix = eigvalues_closeto_zero_matrix(num_dimens);

eigs = eig(matrix);

ratioeig = max(eigs)/min(eigs);

rattable(j,2) = ratioeig;

num_mat = 5; % number of matrices

mat_array = zeros(num_dimens,num_dimens,num_mat);

noise = [0.00001,0,0; % noise matrice 
        0,0,0;
        0,0,0];

% create positive-defined matrices

for i= 1:num_mat
    
    isposdef = 0;

    while isposdef == 0

        matrix = matrix + noise;

        d = eig(matrix);

        isposdef = all(d > 0);

    end

    mat_array(:,:,i) = matrix;

end

T1 = mat_array(:,:,1);

T2 = random_positivedefined_matrices(num_dimens);

T3 = random_positivedefined_matrices(num_dimens);

matrices = zeros(num_dimens,num_dimens,3);

matrices(:,:,1) =  T1;
matrices(:,:,2) =  T2;
matrices(:,:,3) =  T3;


r_m = compute_riepro_mean(matrices,'A',10000);


table = zeros(num_mat-1,3);

meanarray_r = zeros(num_dimens,num_dimens,num_mat-1);

meanarray_r(:,:,1) = r_m;


% compute distance between means

for i = 2:num_mat

    T1 = mat_array(:,:,i);

    matrices(:,:,1) =  T1;


    r_m_t = compute_riepro_mean(matrices,'A',10000);
    
    meanarray_r(:,:,i) = r_m_t;

   

    table(i-1,1)= i-1;

    table(i-1,2)= rie_dis(r_m_t,r_m);


    r_m = r_m_t;


end


    if isreal(meanarray_r) == false
        f = f + 1;
        rattable(j,3) = 1
    else
        rattable(j,3) = 0
    end
end 

complratio = f/tims