
num_dimens = 5;
num_matri = 1000;


sigma1 = random_positivedefined_matrices(num_dimens);

sigma2 = random_positivedefined_matrices(num_dimens);

sigma3 = random_positivedefined_matrices(num_dimens);

dis = rie_dis(sigma2,sigma1);

iter = 100;

for i = 1:iter
  
    dis_1 = rie_dis(sigma3,sigma1);
    
    if dis_1 <dis
         sigma2 = sigma3;
         dis = dis_1;
    end
    
    sigma3 = random_positivedefined_matrices(num_dimens);
    
end

data1 = create_SPD(sigma1,num_dimens,num_matri);
data2 = create_SPD(sigma2,num_dimens,num_matri);

table_a = zeros(9,13);
table_a(:,1) =[1000, 100, 10,1,0.1,0.01,0.001,0.0001,0.00001];

n = 0;
%% 
for j = [1000, 100, 10,1,0.1,0.01,0.001,0.0001,0.00001] 

    n = n+1;
    
    cov1 = data1+noise_matrix_SPD(num_dimens,num_matri)*j;
    cov2 = data2+noise_matrix_SPD(num_dimens,num_matri)*j;

    size1 = size(cov1);
    size2 = size(cov2);

    rat = 0.8; %% proportion of training samples

    cov1train_num = reshape(randsample(size1(3),round(size1(3)*rat)),1,[]);
    cov2train_num = reshape(randsample(size2(3),round(size2(3)*rat)),1,[]);

    cov1train = cov1(:,:,cov1train_num);
    cov2train = cov2(:,:,cov2train_num);

    cov1va_num = setdiff(1:size1(3),cov1train_num);
    cov2va_num = setdiff(1:size2(3),cov2train_num);


    cov1va = cov1(:,:,cov1va_num);
    cov2va = cov2(:,:,cov2va_num);

    cov1ind = size(cov1va);
    cov2ind = size(cov2va);

    cov_va = cat(3,cov1va,cov2va);

    cov1i= ones(cov1ind(3),1);
    cov2i= ones(cov2ind(3),1)*2;


    covind = cat(1,cov1i,cov2i);%% index of va

    %% mean 1
    tic
    class1mean = compute_riepro_mean(cov1train,'A',10000); 
    class2mean = compute_riepro_mean(cov2train,'A',10000); 


    %% validation

    classout = zeros(length(covind),1);

    for i = 1:length(covind)

       point = cov_va(:,:,i);
         dis1 = rie_dis(class1mean,point);
         dis2 = rie_dis(class2mean,point);
       if dis1 > dis2
           classout(i) = 2;
       else 
           classout(i) = 1;
       end 
    end

    lengt = size(covind);
    s_r_accuracyrate = (lengt(1)-sum(abs(classout-covind)))/lengt(1);
    toc1 = toc;
    
    table_a(n,2) = s_r_accuracyrate;
    table_a(n,3) = toc1;
    %% mean 2
    tic
    class1mean = compute_bwpro_mean(cov1train); 
    class2mean = compute_bwpro_mean(cov2train); 


    %% validation

    classout = zeros(length(covind),1);

    for i = 1:length(covind)

       point = cov_va(:,:,i);
         dis1 = compute_W_distance(class1mean,point);
         dis2 = compute_W_distance(class2mean,point);
       if dis1 > dis2
           classout(i) = 2;
       else 
           classout(i) = 1;
       end 
    end

    lengt = size(covind);
    s_w_accuracyrate = (lengt(1)-sum(abs(classout-covind)))/lengt(1);
    toc2 = toc;
    table_a(n,4) = s_w_accuracyrate;
    table_a(n,5) = toc2;
    
    
    %% mean 3
    eps = 1e-2;
    tic
    class1mean = compute_Mean_inductive_rie(cov1train,eps); 
    class2mean = compute_Mean_inductive_rie(cov2train,eps); 


    %% validation

    classout = zeros(length(covind),1);

    for i = 1:length(covind)

       point = cov_va(:,:,i);
         dis1 = rie_dis(class1mean,point);
         dis2 = rie_dis(class2mean,point);
       if dis1 > dis2
           classout(i) = 2;
       else 
           classout(i) = 1;
       end 
    end

    lengt = size(covind);
    i_r_accuracyrate = (lengt(1)-sum(abs(classout-covind)))/lengt(1);
    toc3 = toc;

    table_a(n,6) = i_r_accuracyrate;
    table_a(n,7) = toc3;
    %% mean 4
    tic
    class1mean = compute_Mean_inductive_bw(cov1train,eps);
    class2mean = compute_Mean_inductive_bw(cov1train,eps);


    %% validation

    classout = zeros(length(covind),1);

    for i = 1:length(covind)

       point = cov_va(:,:,i);
         dis1 = compute_W_distance(class1mean,point);
         dis2 = compute_W_distance(class2mean,point);
       if dis1 > dis2
           classout(i) = 2;
       else 
           classout(i) = 1;
       end 
    end

    %%
    lengt = size(covind);
    i_w_accuracyrate = (lengt(1)-sum(abs(classout-covind)))/lengt(1);
    toc4 = toc;
    table_a(n,8) = i_w_accuracyrate;
    table_a(n,9) = toc4;
    
     %% mean 5
    eps = 1e-2;
    tic
    class1mean = compute_riecheap_mean(cov1train); 
    class2mean = compute_riecheap_mean(cov2train); 


    %% validation

    classout = zeros(length(covind),1);

    for i = 1:length(covind)

       point = cov_va(:,:,i);
         dis1 = rie_dis(class1mean,point);
         dis2 = rie_dis(class2mean,point);
       if dis1 > dis2
           classout(i) = 2;
       else 
           classout(i) = 1;
       end 
    end

    lengt = size(covind);
    c_r_accuracyrate = (lengt(1)-sum(abs(classout-covind)))/lengt(1);
    toc5 = toc;
    
    table_a(n,10) = c_r_accuracyrate;
    table_a(n,11) = toc5;

    %% mean 6
    tic
    class1mean = compute_bwcheap_mean(cov1train); 
    class2mean = compute_bwcheap_mean(cov2train); 


    %% validation

    classout = zeros(length(covind),1);

    for i = 1:length(covind)

       point = cov_va(:,:,i);
         dis1 = compute_W_distance(class1mean,point);
         dis2 = compute_W_distance(class2mean,point);
       if dis1 > dis2
           classout(i) = 2;
       else 
           classout(i) = 1;
       end 
    end

    %%
    lengt = size(covind);
    c_w_accuracyrate = (lengt(1)-sum(abs(classout-covind)))/lengt(1);
    toc6 = toc;
    
    table_a(n,12) = c_w_accuracyrate;
    table_a(n,13) = toc6;


    
end

%% 
x = 1:9;
a1 = table_a(:,2); %a r s
b1 = table_a(:,3); %t r s
a2 = table_a(:,4); %a w s
b2 = table_a(:,5); %t w s

a3 = table_a(:,6); %a r i
b3 = table_a(:,7); %t r i
a4 = table_a(:,8); %a w i
b4 = table_a(:,9); %t w i

a5 = table_a(:,10); %a r c
b5 = table_a(:,11); %t r c
a6 = table_a(:,12); %a w c
b6 = table_a(:,13); %t w c

figure(1);  
plot(x,a1,'-*r', 'linewidth', 1.1)
hold on
set(gca,'FontSize',20);
plot(x,a2,'-*b', 'linewidth', 1.1)
xlabel('Intensity of noise(descending)','fontsize',20)  
ylabel('Accuracy','fontsize',20)  
title('Simple Mean')
legend({'rie','bw'},'Location','northwest')

figure(2);  
plot(x,b1,'-*r', 'linewidth', 1.1)
hold on
set(gca,'FontSize',20);
plot(x,b2,'-*b', 'linewidth', 1.1)
xlabel('Intensity of noise(descending)','fontsize',20)  
ylabel('Time','fontsize',20)  
title('Simple Mean')
legend({'rie','bw'},'Location','northwest')

figure(3);  
plot(x,a3,'-*r', 'linewidth', 1.1)
hold on
set(gca,'FontSize',20);
plot(x,a4,'-*b', 'linewidth', 1.1)
xlabel('Intensity of noise(descending)','fontsize',20)  
ylabel('Accuracy','fontsize',20)  
title('Inductive Mean')
legend({'rie','bw'},'Location','northwest')

figure(4);  
plot(x,b3,'-*r', 'linewidth', 1.1)
hold on
set(gca,'FontSize',20);
plot(x,b4,'-*b', 'linewidth', 1.1)
xlabel('Intensity of noise(descending)','fontsize',20)  
ylabel('Time','fontsize',20)  
title('Inductive Mean')
legend({'rie','bw'},'Location','northwest')

figure(5);  
plot(x,a5,'-*r', 'linewidth', 1.1)
hold on
set(gca,'FontSize',20);
plot(x,a6,'-*b', 'linewidth', 1.1)
xlabel('Intensity of noise(descending)','fontsize',20)  
ylabel('Accuracy','fontsize',20)  
title('Cheap Mean')
legend({'rie','bw'},'Location','northwest')

figure(6);  
plot(x,b5,'-*r', 'linewidth', 1.1)
hold on
set(gca,'FontSize',20);
plot(x,b6,'-*b', 'linewidth', 1.1)
xlabel('Intensity of noise(descending)','fontsize',20)  
ylabel('Time','fontsize',20)  
title('Cheap Mean')
legend({'rie','bw'},'Location','northwest')