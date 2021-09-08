

num_dimens = 10; % dimensions of matrices

sigma1 = random_positivedefined_matrices(num_dimens);
eig(sigma1);
df = 2*num_dimens;
[W,D] = wishrnd(sigma1,df);

ran = 10:10:300; % sequence of number of matrices


timetable = zeros(length(ran),6);
iterate = zeros(length(ran),1);
f = 1;

for n = ran
    num_matri = n;
    
    data1 = zeros(num_dimens,num_dimens,num_matri);
    for i= 1:num_matri
        F = wishrnd(sigma1,df,D)/df;
        data1(:,:,i) = F;
    end
    tic
    [mean_w,itera] = compute_Mean1_inductive(data1,eps); 
    tim = toc;
    
    tic
    mean_m = compute_riemannian_mean(data1,"A",15);
    tim_m = toc;
    
    matabs = abs(mean_w-mean_m);
    matmax = max(matabs(:));
    
    r_dis = compute_riemannian_distance(mean_m,mean_w,"A");
    w_dis = compute_W_distance(mean_m,mean_w);
    
    iterate(f) = itera;
    
    timetable(f,1) = n;
    timetable(f,2) = tim;
    timetable(f,3) = tim_m;
    timetable(f,4) = matmax;
    timetable(f,5) = w_dis;
    timetable(f,6) = r_dis;
    f = f+1;
end

x = timetable(:,1);
a = timetable(:,2);
b = timetable(:,3);
c = timetable(:,4);
d = timetable(:,5);
e = timetable(:,6);

figure(1);  
plot(x,a,'-*b', 'linewidth', 1.1)
hold on
plot(x,b,'-or', 'linewidth', 1.1)
set(gca,'FontSize',20);
xlabel('Number of Matrices','fontsize',20)  
ylabel('Running Time','fontsize',20)  
legend({'Bw-distance','R-distance'},'fontsize',15,'Location','northwest');

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
legend({'Bw-distance','R-distance'},'fontsize',15,'Location','northwest');