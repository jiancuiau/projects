a
num_dimens = 3;

C = eigvalues_closeto_zero_matrix(num_dimens);

[V,D] = eig(C);

E = random_positivedefined_matrices(num_dimens)*1e-6;

F = E + C;

%% distance

num_dimens = 3;
E1 = eigvalues_closeto_zero_matrix(num_dimens);
noise = [0.00001,0,0;
        0,0,0;
        0,0,0];

E2 = E1  + noise;

rie_distance1 = rie_dis(E1,E2);

Bw_distance1 = compute_W_distance(E1,E2);

E3 = random_positivedefined_matrices(num_dimens);

noise = [0.001,0,0;
        0,0,0;
        0,0,0];

E4 = E3  + noise;

rie_distance2 = rie_dis(E3,E4);

Bw_distance2 = compute_W_distance(E3,E4);




%%

C(1,1) = C(1,1) +1;
[Y, U] = eig(C);



%% 

num_dimens = 3;

T1 = eigvalues_closeto_zero_matrix(num_dimens);

T2 = random_positivedefined_matrices(num_dimens);

T3 = random_positivedefined_matrices(num_dimens);

matrices = zeros(num_dimens,num_dimens,3);

matrices(:,:,1) =  T1;
matrices(:,:,2) =  T2;
matrices(:,:,3) =  T3;

r_m = rie_mean(matrices,10000);
bw_m = compute_cheap_mean(matrices);

noise = [0.001,0,0;
        0,0,0;
        0,0,0];

iter = 10;

table = zeros(iter,3);

for i = 1:iter
    matrices(:,:,1)= matrices(:,:,1)+noise;
    
    r_m_t = rie_mean(matrices,10000);
    bw_m_t = compute_cheap_mean(matrices);
    
    table(i,1)= i;
    table(i,2)= rie_dis(r_m_t,r_m);
    table(i,3)= compute_W_distance(bw_m_t,bw_m);
    
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
ylabel('distance','fontsize',20)  

figure(2);  
plot(x,b,'-*b', 'linewidth', 1.1)
set(gca,'FontSize',20);
xlabel('Number','fontsize',20)  
ylabel('distance','fontsize',20) 



%% test for path

matrix = eigvalues_closeto_zero_matrix(num_dimens);

%%

num_mat = 100;

mat_array = zeros(num_dimens,num_dimens,num_mat);

noise = [0.00001,0,0;
        0,0,0;
        0,0,0];


for i= 1:100
    
    isposdef = 0;

    while isposdef == 0

        matrix = matrix + noise;

        d = eig(matrix);

        isposdef = all(d > 0);

    end

    mat_array(:,:,i) = matrix;

end