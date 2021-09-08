
ran = 3:20; % sequence of dimensions of matrices
num_matri = 300; % number of matrices
eps = 1e-2; % w-eps

timetable = zeros(length(ran),3);
iterate = zeros(length(ran),1);
f = 1;

for n = ran

    sigma1 = random_positivedefined_matrices(n);
    eig(sigma1);

    df = 2*n;

    [W,D] = wishrnd(sigma1,df);

    data1 = zeros(n,n,num_matri);
    for i= 1:num_matri
        F = wishrnd(sigma1,df,D)/df;
        data1(:,:,i) = F;
    end

    tic
    [mean_w,itera] = compute_Mean_inductive(data1,1000000,eps);
    tim = toc;
    
    tic
    mean_m = compute_riemannian_mean(data1,"A",15);
    tim_m = toc;
    
    iterate(f) = itera;
    
    timetable(f,1) = n;
    timetable(f,2) = tim;
    timetable(f,3) = tim_m;
    f = f+1;
end

x = timetable(:,1);
a = timetable(:,2);
b = timetable(:,3);

plot(x,a,'-*b',x,b,'-or')
title('Number of dimesions VS Time for M and W methods')
xlabel('Number of dimesions')  
ylabel('Time')  
legend({'reiman_w','reiman_m'},'Location','northwest'); 
