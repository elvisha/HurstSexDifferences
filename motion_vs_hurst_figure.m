function motion_vs_hurst_figure(fd_1, fd_2, fd_3, fd_4, fd_start, fd_end, he_1, he_2, he_3, he_4, he_start, he_end);
%generates figure comparing motion parameters between males and females for
%the 4 scans

%fd_1, fd_2, fd_3, fd_4 are each nxm matrices where n is number of subjects
%and m is total time points; each matrix is for a separate scan
%fd_start and fd_end are the start and end time points you want to use

%he_1, he_2, he_3, he_4 are each nxm matrices where n is number of subjects
%and m is total segments; each matrix is for a separate scan
%he_start and he_end are integers that define the start and end segments you 
%want to use to compute he for ech subject


%load in subject sex data
subj_sex=load('subj_sex.txt');

%create matrices for male fd and he data for each scan
male=subj_sex==1;
male_fd1 = cell(size(male));
male_fd1 = fd_1(male,fd_start:fd_end);
male_fd2 = cell(size(male));
male_fd2 = fd_2(male,fd_start:fd_end);
male_fd3 = cell(size(male));
male_fd3 = fd_3(male,fd_start:fd_end);
male_fd4 = cell(size(male));
male_fd4 = fd_4(male,fd_start:fd_end);
male_he1 = cell(size(male));
male_he1 = he_1(male,he_start:he_end);
male_he2 = cell(size(male));
male_he2 = he_2(male,he_start:he_end);
male_he3 = cell(size(male));
male_he3 = he_3(male,he_start:he_end);
male_he4 = cell(size(male));
male_he4 = he_4(male,he_start:he_end);

%create matrices for female fd and he data for each scan
female=subj_sex==0;
female_fd1= cell(size(female));
female_fd1 = fd_1(female,fd_start:fd_end);
female_fd2= cell(size(female));
female_fd2 = fd_2(female,fd_start:fd_end);
female_fd3= cell(size(female));
female_fd3 = fd_3(female,fd_start:fd_end);
female_fd4= cell(size(female));
female_fd4 = fd_4(female,fd_start:fd_end);
female_he1 = cell(size(female));
female_he1 = he_1(female,he_start:he_end);
female_he2 = cell(size(female));
female_he2 = he_2(female,he_start:he_end);
female_he3 = cell(size(female));
female_he3 = he_3(female,he_start:he_end);
female_he4 = cell(size(female));
female_he4 = he_4(female,he_start:he_end);

%generate figure
figure()

%generate subplot to plot scatterplot for motion (x axis) vs hurst (y axis)
%for first scan
%males plotted in blue, females in red
subplot(2,2,1)
scatter(mean(male_fd1,2),mean(male_he1,2),'b')
hold on
scatter(mean(female_fd1,2),mean(female_he1,2),'r')
%set x and y axis labels and limits, title
title('REST1 LR')
xlabel('Motion')
ylabel('BMI')
xlim([0,0.8])
%ylim([0.5,0.7])
%calculate correlation coefficient before fd and he for first scan
%include correlation coefficient r value in legend box
[r1,p1]=corrcoef(mean(fd_1(:,fd_start:fd_end),2),mean(he_1(:,he_start:he_end),2),'rows','complete');
lgd = legend('male', 'female', 'Location', 'southeast');
title(lgd,sprintf('r = %f',r1(1,2)))
grid on
grid minor

%generate subplot to plot scatterplot for motion (x axis) vs hurst (y axis)
%for second scan
%males plotted in blue, females in red
subplot(2,2,2)
scatter(mean(male_fd2,2),mean(male_he2,2),'b')
hold on
scatter(mean(female_fd2,2),mean(female_he2,2),'r')
%set x and y axis labels and limits, title
title('REST1 RL')
xlabel('Motion')
ylabel('Hurst Exponent')
xlim([0,0.8])
ylim([0.5,0.7])
%calculate correlation coefficient before fd and he for second scan
%include correlation coefficient r value in legend box
[r2,p2]=corrcoef(mean(fd_2(:,fd_start:fd_end),2),mean(he_2(:,he_start:he_end),2),'rows','complete');
%male_bmi=corrcoef(mean(male_fd1,2),mean(male_he1,2),'rows','complete')
%female_bmi=corrcoef(mean(female_fd1,2),mean(female_he1,2),'rows','complete')
%male_weight=corrcoef(mean(male_fd1,2),mean(male_he3,2),'rows','complete')
%female_weight=corrcoef(mean(female_fd1,2),mean(female_he3,2),'rows','complete')
lgd = legend('male', 'female', 'Location', 'southeast');
title(lgd,sprintf('r = %f',r2(1,2)))
grid on
grid minor

%generate subplot to plot scatterplot for motion (x axis) vs hurst (y axis)
%for third scan
%males plotted in blue, females in red
subplot(2,2,3)
scatter(mean(male_fd3,2),mean(male_he3,2),'b')
hold on
scatter(mean(female_fd3,2),mean(female_he3,2),'r')
%set x and y axis labels and limits, title
title('REST2 LR')
xlabel('Motion')
ylabel('Weight')
xlim([0,0.8])
ylim([0.5,0.7])
%calculate correlation coefficient before fd and he for third scan
%include correlation coefficient r value in legend box
[r3,p3]=corrcoef(mean(fd_3(:,fd_start:fd_end),2),mean(he_3(:,he_start:he_end),2),'rows','complete');
lgd = legend('male', 'female', 'Location', 'southeast');
title(lgd,sprintf('r = %f',r3(1,2)))
grid on
grid minor

%generate subplot to plot scatterplot for motion (x axis) vs hurst (y axis)
%for fourth scan
%males plotted in blue, females in red
subplot(2,2,4)
scatter(mean(male_fd3,2),mean(male_he3,2),'b')
hold on
scatter(mean(female_fd3,2),mean(female_he3,2),'r')
title('REST2 RL')
%set x and y axis labels and limits, title
xlabel('Motion')
ylabel('Hurst Exponent')
xlim([0,0.8])
ylim([0.5,0.7])
%calculate correlation coefficient before fd and he for fourth scan
%include correlation coefficient r value in legend box
[r4,p4]=corrcoef(mean(fd_4(:,fd_start:fd_end),2),mean(he_4(:,he_start:he_end),2),'rows','complete');
lgd = legend('male', 'female', 'Location', 'southeast');
title(lgd,sprintf('r = %f',r4(1,2)))
grid on
grid minor
end
