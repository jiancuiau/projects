num_dimens = 10; % dimensions of matrices

ran = 10:10:300; % sequence of number of matrices

timetable = zeros(length(ran),6);
f = 1;
esp = 1e-6;
itertable = zeros(length(ran),6);

for n = ran
    n
    num_matri = n;
    
    data1 = random_positivedefined_matrices_number(num_dimens,num_matri);
    
    tic
    [mean_1,it1] = compute_bwcheap_mean(data1); 
    tim = toc;
    
    tic
    [mean_2,it2] = compute_riecheap_mean(data1);
    tim_m = toc;
    
    tic
    [mean_3,it3] = compute_Mean_inductive_bw_n(data1,esp);
    tim_m1 = toc;
    
    tic
    [mean_4,it4] = compute_Mean_inductive_rie_n(data1,esp);
    tim_m2 = toc;
    
    tic
    [mean_5,it5] = compute_bwpro_mean(data1);
    tim_m3 = toc;
    
    tic
    [mean_6,it6] = compute_riepro_mean(data1,"A",1000);
    tim_m4 = toc;
    
    
    matabs12 = abs(mean_1-mean_2);
    matmax12 = max(matabs12(:));
    
    dis_2 = compute_riemannian_distance(mean_1,mean_2,"A");
    dis_1 = compute_W_distance(mean_1,mean_2);
    
    matabs34 = abs(mean_3-mean_4);
    matmax34 = max(matabs34(:));
    
    dis_2_34 = compute_riemannian_distance(mean_3,mean_4,"A");
    dis_1_34 = compute_W_distance(mean_3,mean_4);
    
    matabs56 = abs(mean_5-mean_6);
    matmax56 = max(matabs56(:));
    
    dis_2_56 = compute_riemannian_distance(mean_5,mean_6,"A");
    dis_1_56 = compute_W_distance(mean_5,mean_6);
    
    timetable(f,1) = n;
    timetable(f,2) = tim;
    timetable(f,3) = tim_m;
    timetable(f,4) = matmax12;
    timetable(f,5) = dis_1;
    timetable(f,6) = dis_2;
    timetable(f,7) = tim_m1; %indu_bw
    timetable(f,8) = tim_m2; %indu_rie
    timetable(f,9) = matmax34;
    timetable(f,10) = dis_1_34;
    timetable(f,11) = dis_2_34;
    
    timetable(f,12) = tim_m3; %simple projection_bw
    timetable(f,13) = tim_m4; %simple projection_rie
    timetable(f,14) = matmax56;
    timetable(f,15) = dis_1_56;
    timetable(f,16) = dis_2_56;
    
    
    
    itertable(f,1) = it1;
    itertable(f,2) = it2;
    itertable(f,3) = it3;
    itertable(f,4) = it4;
    itertable(f,5) = it5;
    itertable(f,6) = it6;
    
    f = f+1;
end


x = timetable(:,1);
a = timetable(:,2);
b = timetable(:,3);
c = timetable(:,4);
d = timetable(:,5);
e = timetable(:,6);
r = a./b;

a1 = timetable(:,7);
b1 = timetable(:,8);
c1 = timetable(:,9);
d1 = timetable(:,10);
e1 = timetable(:,11);
r1 = a1./b1;
a2 = timetable(:,12);
b2 = timetable(:,13);
c2 = timetable(:,14);
d2 = timetable(:,15);
e2 = timetable(:,16);
r2 = a2./b2;



% cheap graphs
figure(1);  
yyaxis left
plot(x,a,'-*b', 'linewidth', 1.1)
hold on
plot(x,b,'-or', 'linewidth', 1.1)
set(gca,'FontSize',20);
ylabel('Running Time','fontsize',20) 
yyaxis right
plot(x,r,'-g', 'linewidth', 1.1);
set(gca,'FontSize',20);
xlabel('Number of Matrices','fontsize',20)  
ylabel('Ratio of Running Time','fontsize',20)  
legend({'Bw-cheap','Rie-cheap','Ratio of Bw/Rie'},'fontsize',15,'Location','northwest');

figure(2);  
plot(x,c,'-*b', 'linewidth', 1.1)
set(gca,'FontSize',20);
xlabel('Number of Matrices','fontsize',20)  
ylabel('Matrices Difference','fontsize',20) 

figure(3);
plot(x,d,'-*b', 'linewidth', 1.1)
hold on
plot(x,e,'-or', 'linewidth', 1.1)
set(gca,'FontSize',20);
xlabel('Number of Matrices','fontsize',20)  
ylabel('Distance','fontsize',20)  
legend({'Bw-cheap','Rie-cheap'},'fontsize',15,'Location','northwest');

% inductive graphs
figure(4); 
yyaxis left
plot(x,a1,'-*b', 'linewidth', 1.1)
hold on
plot(x,b1,'-or', 'linewidth', 1.1)
set(gca,'FontSize',20);
ylabel('Running Time','fontsize',20) 
yyaxis right
plot(x,r1,'-g', 'linewidth', 1.1);
set(gca,'FontSize',20);
xlabel('Number of Matrices','fontsize',20)  
ylabel('Ratio of Running Time','fontsize',20) 
legend({'Bw-inductive','Rie-inductive','Ratio of Bw/Rie'},'fontsize',15,'Location','northwest');

figure(5);  
plot(x,c1,'-*b', 'linewidth', 1.1)
set(gca,'FontSize',20);
xlabel('Number of Matrices','fontsize',20)  
ylabel('Matrices Difference','fontsize',20) 

figure(6);
plot(x,d1,'-*b', 'linewidth', 1.1)
hold on
plot(x,e1,'-or', 'linewidth', 1.1)
set(gca,'FontSize',20);
xlabel('Number of Matrices','fontsize',20)  
ylabel('Distance','fontsize',20)  
legend({'Bw-inductive','Rie-inductive'},'fontsize',15,'Location','northwest');


% simple projection graphs
figure(7); 
yyaxis left
plot(x,a2,'-*b', 'linewidth', 1.1)
hold on
plot(x,b2,'-or', 'linewidth', 1.1)
set(gca,'FontSize',20);
ylabel('Running Time','fontsize',20) 
yyaxis right
plot(x,r2,'-g', 'linewidth', 1.1);
set(gca,'FontSize',20);
xlabel('Number of Matrices','fontsize',20)  
ylabel('Ratio of Running Time','fontsize',20) 
legend({'Bw-simple','Rie-simple','Ratio of Bw/Rie'},'fontsize',15,'Location','northwest');

figure(8);  
plot(x,c2,'-*b', 'linewidth', 1.1)
set(gca,'FontSize',20);
xlabel('Number of Matrices','fontsize',20)  
ylabel('Matrices Difference','fontsize',20) 

figure(9);
plot(x,d2,'-*b', 'linewidth', 1.1)
hold on
plot(x,e2,'-or', 'linewidth', 1.1)
set(gca,'FontSize',20);
xlabel('Number of Matrices','fontsize',20)  
ylabel('Distance','fontsize',20)  
legend({'Bw-simple','Rie-simple'},'fontsize',15,'Location','northwest');


figure(10);  
plot(itertable,'linewidth', 1.1)
set(gca,'FontSize',20);
xlabel('Number of Matrices','fontsize',20)  
ylabel('Number of Iteration','fontsize',20)
legend({'Bw-cheap','Rie-cheap','Bw-indu','Rie-indu','Bw-pro','Rie-pro'},'fontsize',15,'Location','northeast');