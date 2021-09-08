% graph

x = 1:9;
a1 = table(:,2); %a r s
b1 = table(:,3); %t r s
a2 = table(:,4); %a w s
b2 = table(:,5); %t w s

a3 = table(:,6); %a r i
b3 = table(:,7); %t r i
a4 = table(:,8); %a w i
b4 = table(:,9); %t w i

a5 = table(:,10); %a r c
b5 = table(:,11); %t r c
a6 = table(:,12); %a w c
b6 = table(:,13); %t w c





figure(1);  
plot(x,a1,'-*b', 'linewidth', 1.1)
set(gca,'FontSize',20);
plot(x,a2,'-*b', 'linewidth', 1.1)
xlabel('Intensity of noise(descending)','fontsize',20)  
ylabel('Accuracy','fontsize',20)  
title('Simple Mean')
legend({'rie','bw'},'Location','northwest')

figure(2);  
plot(x,b1,'-*b', 'linewidth', 1.1)
set(gca,'FontSize',20);
plot(x,b2,'-*b', 'linewidth', 1.1)
xlabel('Intensity of noise(descending)','fontsize',20)  
ylabel('Time','fontsize',20)  
title('Simple Mean')
legend({'rie','bw'},'Location','northwest')

figure(3);  
plot(x,a3,'-*b', 'linewidth', 1.1)
set(gca,'FontSize',20);
plot(x,a4,'-*b', 'linewidth', 1.1)
xlabel('Intensity of noise(descending)','fontsize',20)  
ylabel('Accuracy','fontsize',20)  
title('Inductive Mean')
legend({'rie','bw'},'Location','northwest')

figure(4);  
plot(x,b3,'-*b', 'linewidth', 1.1)
set(gca,'FontSize',20);
plot(x,b4,'-*b', 'linewidth', 1.1)
xlabel('Intensity of noise(descending)','fontsize',20)  
ylabel('Time','fontsize',20)  
title('Inductive Mean')
legend({'rie','bw'},'Location','northwest')

figure(5);  
plot(x,a5,'-*b', 'linewidth', 1.1)
set(gca,'FontSize',20);
plot(x,a6,'-*b', 'linewidth', 1.1)
xlabel('Intensity of noise(descending)','fontsize',20)  
ylabel('Accuracy','fontsize',20)  
title('Cheap Mean')
legend({'rie','bw'},'Location','northwest')

figure(6);  
plot(x,b5,'-*b', 'linewidth', 1.1)
set(gca,'FontSize',20);
plot(x,b6,'-*b', 'linewidth', 1.1)
xlabel('Intensity of noise(descending)','fontsize',20)  
ylabel('Time','fontsize',20)  
title('Cheap Mean')
legend({'rie','bw'},'Location','northwest')



%% sp

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
