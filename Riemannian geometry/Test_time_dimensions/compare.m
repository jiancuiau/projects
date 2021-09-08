%%
normtableall = zeros(4,4,1);

normtable = zeros(4,4);

itertable = zeros(4,6);

seq = [15,25,35,45];
for j =1:1
for i = 1:4

num_dimens = seq(i);
seq(i)
num_matri = 10;

matrices = random_positivedefined_matrices_number(num_dimens,num_matri);

tic
[mean1,it1] = compute_bwcheap_mean(matrices);
[mean2,it2] = compute_riecheap_mean(matrices);
[mean3,it3] = compute_bwpro_mean(matrices);
[mean4,it4] = compute_riepro_mean(matrices,'A',1000);
[mean5,it5] = compute_Mean_inductive_bw_n(matrices,1e-6);
[mean6,it6] = compute_Mean_inductive_rie_n(matrices,1e-6);
toc


bwcheapnorm = norm(mean1-mean5)/norm(mean5);
riecheapnorm = norm(mean2-mean6)/norm(mean6);
bwpronorm = norm(mean3-mean5)/norm(mean5);
riepronorm = norm(mean4-mean6)/norm(mean6);

normtable(i,1) = bwcheapnorm;
normtable(i,2) = riecheapnorm;
normtable(i,3) = bwpronorm;
normtable(i,4) = riepronorm;

    for h = 1:6
        itertable(i,h) = eval(['it',num2str(h)]);
    end
end
normtableall(:,:,j) = normtable;

end
%%
seq = [3,5,10,20,30,40];


normre = normre(1:6,:);

a = normre(:,1) ;
b = normre(:,2) ;
c = normre(:,3) ;
d = normre(:,4) ;

%%
figure(1);  
plot(seq,a,'-*b', 'linewidth', 1.1)
hold on
plot(seq,b,'-or', 'linewidth', 1.1)
set(gca,'FontSize',20);
xlabel('Number of Dimensions','fontsize',20)  
ylabel('Norm','fontsize',20)  
legend({'Bw-cheap','Rie-cheap'},'fontsize',15,'Location','northwest');



figure(2);  
plot(seq,c,'-*b', 'linewidth', 1.1)
hold on
plot(seq,d,'-or', 'linewidth', 1.1)
set(gca,'FontSize',20);
xlabel('Number of Dimensions','fontsize',20)  
ylabel('Norm','fontsize',20)  
legend({'Bw-pro','Rie-pro'},'fontsize',15,'Location','northwest');

%%


num_dimens = 5;
num_matri = 10;

matrices = random_positivedefined_matrices_number(num_dimens,num_matri);
%%

tic
mean1 = compute_bwcheap_mean(matrices);
mean2 = compute_riecheap_mean(matrices);
mean3 = compute_bwpro_mean(matrices);
mean4 = compute_riepro_mean(matrices,'A',1000);
mean5 = compute_Mean_inductive_bw_n(matrices,1e-2);
mean6 = compute_Mean_inductive_rie_n(matrices,1e-2);
mean7 = compute_Mean_inductive_bw_n(matrices,1e-6);
mean8 = compute_Mean_inductive_rie_n(matrices,1e-6);
toc

%%

bwcheapnorm = norm(mean1-mean5)/norm(mean5);
riecheapnorm = norm(mean2-mean6)/norm(mean6);
bwpronorm = norm(mean3-mean5)/norm(mean5);
riepronorm = norm(mean4-mean6)/norm(mean6);

bwcheapnorm6 = norm(mean1-mean7)/norm(mean7);
riecheapnorm6 = norm(mean2-mean8)/norm(mean8);
bwpronorm6 = norm(mean3-mean7)/norm(mean7);
riepronorm6 = norm(mean4-mean8)/norm(mean8);

%%
meanstore = zeros(5,5,8);

meanstore(:,:,1)= mean1;
meanstore(:,:,2)= mean2;
meanstore(:,:,3)= mean3;
meanstore(:,:,4)= mean4;
meanstore(:,:,5)= mean5;
meanstore(:,:,6)= mean6;
meanstore(:,:,7)= mean7;
meanstore(:,:,8)= mean8;

%%
csvwrite('1.csv',mean1);
csvwrite('2.csv',mean2);
csvwrite('3.csv',mean3);
csvwrite('4.csv',mean4);
csvwrite('5.csv',mean5);
csvwrite('6.csv',mean6);

%% 

csvwrite('7.csv',meanstore);




%%



%%
normtableall = zeros(4,4,1);

normtable = zeros(4,4);

itertable = zeros(4,6);

seq = [35];
for j =1:1
    for i = 1:1

    num_dimens = seq(i);
    seq(i)
    num_matri = 3;

    matrices = random_positivedefined_matrices_number(num_dimens,num_matri);

    tic
    
    spd_matrices = matrices;
    
    max_iter = 1000;
    
    [dims,~,num_spd] = size(spd_matrices);
    
    % compute inital mean
    M  = mean(spd_matrices, 3);
    n = 0;      
    meanarray = zeros(num_dimens,num_dimens,20);
            for ite_th = 1 : max_iter
                n = n + 1;
                A = M ^ (1/2);      %-- A = C^(1/2)
                B = A ^ (-1);       %-- B = C^(-1/2)

                S = zeros(size(M));
                for j_th = 1 : num_spd
                    C = spd_matrices(:,:,j_th);
                    S = S + A * logm(B * C * B) * A;
                end
                S = S / num_spd;

                M = A * expm(B * S * B) * A; 
                
                M = (M + M.')/2;
                
                meanarray(:,:,n) = M;
                
                eps = norm(S, 'fro');
                if (eps < 1e-6)
                    break;
                end
            end
            
        
    dis = eps;
    mean_center = M;

    toc



    end
normtableall(:,:,j) = normtable;

end

distable = zeros(5,1);
p=0;

for i = 1+p:ite_th-1+p
    
dis = rie_dis(meanarray(:,:,i),meanarray(:,:,i+1));

distable(i) = dis;

end

plot(distable)


%%
save('iterall.mat')


%%
num_dimens = 25;
    num_matri = 10;

    matrices = random_positivedefined_matrices_number(num_dimens,num_matri);
%%
    tic

    [mean4,it4] = compute_riepro_mean(matrices,'A',2000);

    toc
    
%%

distable = zeros(5,1);
p=0;

for i = 1+p:ite_th-1+p
    
dis = rie_dis(meanarray(:,:,i),meanarray(:,:,i+1));

distable(i) = dis;

end

plot(distable)

%%
meanarray(:,:,402)
meanarray(:,:,403)
meanarray(:,:,404)



%%
normtableall = zeros(4,4,1);

normtable = zeros(4,4);

itertable = zeros(4,6);

seq = [3,5,8,10,15,20,25];

seqn = [3,5,10,15,20,30,35];

cyctab =zeros(7,7);
for h = 1:7
    cycvec = zeros(7,1);
    h
    for i =1:7

        cyc = 0;
        for j = 1:100

        num_dimens = seq(i);
        seq(i);
        num_matri = seqn(h);

        matrices = random_positivedefined_matrices_number(num_dimens,num_matri);

        spd_matrices = matrices;

        max_iter = 1000;

        [dims,~,num_spd] = size(spd_matrices);

        % compute inital mean
        M  = mean(spd_matrices, 3);
        n = 0;      
        meanarray = zeros(num_dimens,num_dimens,20);
                for ite_th = 1 : max_iter
                    n = n + 1;
                    A = M ^ (1/2);      %-- A = C^(1/2)
                    B = A ^ (-1);       %-- B = C^(-1/2)

                    S = zeros(size(M));
                    for j_th = 1 : num_spd
                        C = spd_matrices(:,:,j_th);
                        S = S + A * logm(B * C * B) * A;
                    end
                    S = S / num_spd;

                    M = A * expm(B * S * B) * A; 

                    M = (M + M.')/2;

                    meanarray(:,:,n) = M;

                    eps = norm(S, 'fro');
                    if (eps < 1e-6)
                        break;
                    end
                end


        dis = eps;

        if dis > 0.1
            cyc = cyc + 1;
        end
        mean_center = M;

        end

        cycvec(i) = cyc;
    end
 cyctab(:,h) = cycvec;
end


%%
distable = zeros(5,1);
p=0;

for i = 1+p:ite_th-1+p
    
dis = rie_dis(meanarray(:,:,i),meanarray(:,:,i+1));

distable(i) = dis;

end

plot(distable)




%% inductive - simple projection verify

max_iter =1000;

seqd = [3,5,8,10,15,20,25,40,60,80,100];

seqn = [3,5,10,15,20,30,35,40,60,80,100];


itertable = zeros(length(seqn),length(seqd));

distable = zeros(length(seqn),length(seqd));

for j = 1:length(seqd)
j
num_dimens = seqd(j);

itevec = zeros(length(seqn),1);

disvec = zeros(length(seqn),1);

    for i = 1:length(seqn)


    num_matri = seqn(i);

    es = 1e-2;

    matrices = random_positivedefined_matrices_number(num_dimens,num_matri);

    [mean5,it5] = compute_Mean_inductive_rie_n(matrices,es);


        [dims,~,num_spd] = size(matrices);

        % compute inital mean
        M  = mean5;
        n = 0;

    for ite_th = 1 : max_iter

        n = n + 1;
        A = M ^ (1/2);      %-- A = C^(1/2)
        B = A ^ (-1);       %-- B = C^(-1/2)

        S = zeros(size(M));

        for j_th = 1 : num_spd
            C = matrices(:,:,j_th);
            S = S + A * logm(B * C * B) * A;
        end
        S = S / num_spd;

        M = A * expm(B * S * B) * A; 

        M = (M + M.')/2;

        eps = norm(S, 'fro');
        if (eps < es)
            break;
        end

     end

    dis = eps;

    mean_center = M;

    dismean = rie_dis(mean5,M);
    
    itevec(i,1) = n;
    disvec(i,1) = dismean;
    
    end
    
    
    itertable(:,j) = itevec;
    distable(:,j) = disvec;
    
end

%% inductive - BW simple projection

max_iter =1000;

seqd = [3,5,8,10,15,20,25,40,60,80,100];

seqn = [3,5,10,15,20,30,35,40,60,80,100];


itertable = zeros(length(seqn),length(seqd));

distable = zeros(length(seqn),length(seqd));

for j = 1:length(seqd)
    
j

num_dimens = seqd(j);

itevec = zeros(length(seqn),1);

disvec = zeros(length(seqn),1);

    for i = 1:length(seqn)

    num_matri = seqn(i);

    es = 1e-4;

    matrices = random_positivedefined_matrices_number(num_dimens,num_matri);

    [mean5,it5] = compute_Mean_inductive_bw_n(matrices,es);


 % compute inital mean
    M = mean5;
    
    [~,num_de,num_mat] = size(matrices);
    
% precision of algorithm
    dis = 10;
    
    n = 0;
    
    while dis > es

        n = n + 1;
        M_o = M;
        [V,D] = eig(M);
        V_1 = V^(-1);
        eigv = eig(M);
        s_matrix = zeros(num_de,num_de,num_mat);
        
        % projection to tangent space
        for h = 1:1:num_mat
            
            p = matrices(:,:,h);
            s_matrix(:,:,h) = (M * p)^(1/2)+(p * M)^(1/2) - 2 * M;
            
        end
        
        S = mean(s_matrix,3);

        W = zeros(num_de,num_de);

        for l = 1:1:num_de
            
            for r = 1:1:num_de
                va =(eigv(l)+eigv(r))^(-1);
                W(l,r) = va;
            end
            
        end
        
        % projection to manifold
        
        M = V*(W.*(V_1*S*V+2*D))*D*(W.*(V_1*S*V+2*D))*V_1;

        dis = compute_W_distance(M,M_o);
    
     end


    mean_center = M;

    dismean = compute_W_distance(mean5,M);
    
    itevec(i,1) = n;
    
    disvec(i,1) = dismean;
    
    end
    
    
    itertable(:,j) = itevec;
    
    distable(:,j) = disvec;
    
end

%% BW in-4 sp-2



max_iter =1000;

seqd = [3,5,8,10,15,20,25,40,60,80,100];

seqn = [3,5,10,15,20,30,35,40,60,80,100];


itertable = zeros(length(seqn),length(seqd));

distable = zeros(length(seqn),length(seqd));

for j = 1:length(seqd)
    
j

num_dimens = seqd(j);

itevec = zeros(length(seqn),1);

disvec = zeros(length(seqn),1);

    for i = 1:length(seqn)

    num_matri = seqn(i);

    es_in = 1e-4;
    es_sp = 1e-2;

    matrices = random_positivedefined_matrices_number(num_dimens,num_matri);

    [mean5,it5] = compute_Mean_inductive_bw_n(matrices,es_in);


 % compute inital mean
    M = mean5;
    
    [~,num_de,num_mat] = size(matrices);
    
% precision of algorithm
    dis = 10;
    
    n = 0;
    
    while dis > es_sp

        n = n + 1;
        M_o = M;
        [V,D] = eig(M);
        V_1 = V^(-1);
        eigv = eig(M);
        s_matrix = zeros(num_de,num_de,num_mat);
        
        % projection to tangent space
        for h = 1:1:num_mat
            
            p = matrices(:,:,h);
            s_matrix(:,:,h) = (M * p)^(1/2)+(p * M)^(1/2) - 2 * M;
            
        end
        
        S = mean(s_matrix,3);

        W = zeros(num_de,num_de);

        for l = 1:1:num_de
            
            for r = 1:1:num_de
                va =(eigv(l)+eigv(r))^(-1);
                W(l,r) = va;
            end
            
        end
        
        % projection to manifold
        
        M = V*(W.*(V_1*S*V+2*D))*D*(W.*(V_1*S*V+2*D))*V_1;

        dis = compute_W_distance(M,M_o);
    
     end


    mean_center = M;

    dismean = compute_W_distance(mean5,M);
    
    itevec(i,1) = n;
    
    disvec(i,1) = dismean;
    
    end
    
    
    itertable(:,j) = itevec;
    
    distable(:,j) = disvec;
    
end


%% Rie in-4 sp-2

max_iter =1000;

seqd = [3,5,8,10,15,20,25,40,60,80,100];

seqn = [3,5,10,15,20,30,35,40,60,80,100];


itertable = zeros(length(seqn),length(seqd));

distable = zeros(length(seqn),length(seqd));

for j = 1:length(seqd)
j
num_dimens = seqd(j);

itevec = zeros(length(seqn),1);

disvec = zeros(length(seqn),1);

    for i = 1:length(seqn)


    num_matri = seqn(i);

    es_in = 1e-4;
    es_sp = 1e-2;

    matrices = random_positivedefined_matrices_number(num_dimens,num_matri);

    [mean5,it5] = compute_Mean_inductive_rie_n(matrices,es_in);


        [dims,~,num_spd] = size(matrices);

        % compute inital mean
        M  = mean5;
        n = 0;

    for ite_th = 1 : max_iter

        M_0 = M;
        n = n + 1;
        A = M ^ (1/2);      %-- A = C^(1/2)
        B = A ^ (-1);       %-- B = C^(-1/2)

        S = zeros(size(M));

        for j_th = 1 : num_spd
            C = matrices(:,:,j_th);
            S = S + A * logm(B * C * B) * A;
        end
        S = S / num_spd;

        M = A * expm(B * S * B) * A; 

        M = (M + M.')/2;

        eps = rie_dis(M,M_0);
        if (eps < es_sp)
            break;
        end

     end

    dis = eps;

    mean_center = M;

    dismean = rie_dis(mean5,M);
    
    itevec(i,1) = n;
    disvec(i,1) = dismean;
    
    end
    
    
    itertable(:,j) = itevec;
    distable(:,j) = disvec;
    
end