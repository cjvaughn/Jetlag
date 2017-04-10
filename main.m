%Jetlag
%Note: t is in hours

% clear memory
clearvars
% start the timer
tic
% where to save the data
jobstring='april_7_Ex1'

%% Initializing Parameters:
% for checking if sum is 1, and alpha<alpha_max, V>0
threshold=10^(-5)

omega_0=2*pi/24.5
omega_S=2*pi/24
p=0 %time zone shift %ToDo: make function of t
kappa=1
R=1 %1/(2*omega_S^4) %coupling cost strength %ToDo: set below R_c(0)
F=0 %sun cost strength

% number of times to iterate between HJB and Kolmogorov
num_iterations=40 %ToDo

% this is the maximum alpha that can be reached to meet
% the stability condition
alpha_max=1
alpha_min=-alpha_max

% dynamics: dtheta=(omega_0+alpha)*dt+sigma*dW
sigma=0.1

% number of time steps
num_time_points=1001
% number of grid points in theta
num_x=128;
delta_x=2*pi/num_x;
x_min=0
x_max=2*pi-delta_x;

% grid size in time
delta_t=0.5*1/(sigma^2/(delta_x)^2+((omega_0+alpha_max)/(delta_x)));
delta_t
% finite time horizon, T
T=(num_time_points-1)*delta_t;
T

t_grid=linspace(0,T,num_time_points);
x_grid=linspace(x_min,x_max,num_x);

%% Initializing iterating:

% counter for iterations between HJB and Kolmogorov
K=1;
% initialize mu
mu=zeros(num_time_points,num_x);
%uniform
%mu(:,:)=1/(num_x*delta_x);
%sine
% for i=1:num_time_points
%     for j=1:num_x
%         mu(i,j)=(1+sin(omega_S*t_grid(i)-x_grid(j)))/(num_x*delta_x);
%     end
% end
%cosine
% for i=1:num_time_points
%     for j=1:num_x
%         mu(i,j)=(1+cos(omega_S*t_grid(i)-x_grid(j)))/(num_x*delta_x);
%     end
% end
%dirac at 0
% mu(:,1)=1/delta_x;
%dirac at omega_0*t
for i=1:num_time_points
    index=round(omega_0*t_grid(i)/delta_x)+1;
    mu(i,index)=1/delta_x;
end
value=sum(sum(mu(1,:)))*delta_x;
if value>1+threshold || value<1-threshold
    'Oh no!!!!! Initial mu does not sum to 1'
    value
end

% used to determine when we have converged
mu_diff_max=1;
% alpha from previous iteration
old_alpha=zeros(num_time_points,num_x);

%% Iteration:
while(mu_diff_max>0 && K<num_iterations+1) %K<2 means 1 iteration
old_mu=mu;
'Part B: HJB'
%% Given mu, solve for V (explicitly backwards in time)

% max_alpha is the largest alpha calculated
% (which is different from alpha_max, which is the largest allowed alpha)
max_alpha=0;
V=zeros(num_time_points,num_x);
V(num_time_points,:)=0;

max_diff_alpha=0;
for counter=1:num_time_points-1
    n=num_time_points-counter;
    t_n=t_grid(n);
    
    mu_curr=squeeze(mu(n,:));
    V_curr=squeeze(V(n+1,:));

    V_xx=shift(V_curr,1)-2*V_curr+shift(V_curr,-1); %centered

    right_V_x=(shift(V_curr,1)-V_curr);
    left_V_x=(V_curr-shift(V_curr,-1));
    right=omega_0-1/R*right_V_x/delta_x;
    left=omega_0-1/R*left_V_x/delta_x;
    V_x=zeros(1,num_x);
    V_x(left>0 & right>0)=right_V_x(left>0 & right>0);
    V_x(left<0 & right<0)=left_V_x(left<0 & right<0);
    
    alpha=-V_x/(R*delta_x);
    
    c_bar=zeros(1,num_x); %ToDo: use convolution to be faster
    for i=1:num_x
        for j=1:num_x
            c_bar(i)=c_bar(i)+1/2*sin((x_grid(j)-x_grid(i))/2)^2*mu_curr(j)*delta_x;
        end
    end
    c_sun=1/2*sin(omega_S*t_n+p-x_grid).^2;

    if K>1
        V(n,:)=V_curr+delta_t*(sigma^2/2*V_xx/(delta_x)^2+(omega_0+alpha).*V_x/delta_x+R/2*(squeeze(old_alpha(n,:)).*alpha)+kappa*c_bar+F*c_sun);
    else
        V(n,:)=V_curr+delta_t*(sigma^2/2*V_xx/(delta_x)^2+(omega_0+alpha).*V_x/delta_x+R/2*(alpha.*alpha)+kappa*c_bar+F*c_sun);
    end
    
    if K>1
        alpha_diff=abs(squeeze(old_alpha(n,:))-alpha);
        alpha_diff_real=zeros(num_x,1);
        alpha_diff_real(mu_curr>0)=alpha_diff(mu_curr>0);
        max_curr=max(alpha_diff_real);
        max_diff_alpha=max(max_curr,max_diff_alpha);
    end

    max_alpha=max(max_alpha,max(abs(alpha)));
    old_alpha(n,:)=alpha;
end
%Checking if the solution is valid:
if K>1
    'Difference in alpha'
    max_diff_alpha
end
'Part B Max Alpha'
max_alpha
if max_alpha>alpha_max+threshold
    'Oh no!!!!! Part B max_alpha>alpha_max'
    max_alpha
end
'Minimum of V'
value=min(min(V(:,:)))
if value<-threshold
    'Oh no!!!!! Negatives in V'
    value
end
'Maximum of V'
value=max(max(V(:,:)))

'Part C: Kolmogorov'
%% Given V, solve for mu (explicitly forward in time)
max_alpha=0;
mu=zeros(num_time_points,num_x);
%uniform
%mu(1,:)=1/(num_x*delta_x);
%sine
%mu(1,:)=(1+sin(x_grid))/(num_x*delta_x);
%cosine
% mu(1,:)=(1+cos(x_grid))/(num_x*delta_x);
%dirac at 0
% mu(1,1)=1/delta_x;
%dirac at omega_0*t
mu(1,1)=1/delta_x;

for n=1:num_time_points-1
    t_n=t_grid(n);
        
    mu_curr=squeeze(mu(n,:));
    V_curr=squeeze(V(n,:));
    
    mu_xx=shift(mu_curr,1)-2*mu_curr+shift(mu_curr,-1); %centered

    alpha=squeeze(old_alpha(n,:));
    drift=omega_0+alpha;
    
    drift_minus=zeros(1,num_x);
    drift_plus=zeros(1,num_x);
    drift_minus(drift<0)=-drift(drift<0);
    drift_plus(drift>0)=drift(drift>0);
    
    term_a=shift(drift_minus,1);
    term_b=shift(drift_plus,-1);
    term_c=-(drift_minus+drift_plus);
    
    mu(n+1,:)=mu_curr+delta_t*(sigma^2/2*mu_xx/(delta_x)^2+term_a.*shift(mu_curr,1)/delta_x+term_b.*shift(mu_curr,-1)/delta_x+term_c.*mu_curr/delta_x);

    max_alpha=max(max_alpha,max(abs(alpha)));
end

%Checking if the solution is valid:
integral_values=zeros(num_time_points,1);
for n=1:num_time_points
    integral=sum(mu(n,:))*delta_x;
    integral_values(n)=integral;
    if integral>1+threshold || integral<1-threshold
        'Oh no!!!!! Sum is not 1!!!!'
        integral
    end
end
'Part C Max Alpha'
max_alpha
if max_alpha>alpha_max+threshold
    'Oh no!!!!! Part C max_alpha>alpha_max'
    max_alpha
end
%'Minimum Density in Part C'
value=min(min(mu(:,:)));
if value<-threshold
    'Oh no!!!!! Negatives in mu'
    value
end
mu_diff=abs(mu-old_mu);
mu_diff_max=max(max(mu_diff(:,:)))*delta_x;
value=mu_diff_max;
K
'Difference in mu'
value
mu_diff_frac=mu_diff./abs(mu);
mu_diff_frac(abs(mu)<10^(-10) & abs(old_mu)<10^(-10))=-1;
value=max(max(mu_diff_frac(:,:)));
'Largest Fractional Difference'
value
K=K+1;

%% Final calculations and save data







end %This ends the while loop

%% Final calculations and save data
% stop the timer
timer=toc
'Done'












