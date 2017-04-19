num_x=256;
delta_x=2*pi/num_x;
x_min=0
x_max=2*pi-delta_x;
x_grid=linspace(x_min,x_max,num_x);

thing=V_hourly;

time=1;
f=thing(time,:);
num_time_points=size(thing(:,1),1);

plot(x_grid,f)
title('mu(t,theta), R=1, F=1')
xlabel('theta')
ylabel('mu(T,theta)')

jump=1;

while true
    [~,~,b] = ginput(1);
    if b==29
        time=time+jump;
    elseif b==28
        time=time-jump;
    elseif b==27
        break
    end
    if time>0 && time<=num_time_points
        time
        f=thing(time,:);
        
        plot(x_grid,f)
        title('mu(t,theta), R=1, F=1')
        xlabel('theta')
        ylabel('mu(T,theta)')
    end
end