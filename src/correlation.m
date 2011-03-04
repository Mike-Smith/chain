function cor=correlation(s)

a=importdata(s);
b=size(a);
b=b(1);
a=a(:,1:2:3);
dt=a(2,1);
%plot(a(1024*15:1024*15+128,1),a(1024*15:1024*15+128,2));
T=512;
for i=1:T
    if mod(i,1024)==0
        i
    end
    c(i)=0;
    for j=1:b-T
        c(i)=c(i)+a(j,2)*a(j+i,2);
    end
    c(i)=c(i)/(b-T);
end
time=dt:dt:T*dt;
plot(time(1:T),c(1:T),'Color','k')
hold on
b=importdata('b.dat');
plot(b(1:T,1),b(1:T,2),'Color','r');
plot(time(1:64:T),(8.314472477042e-3)*(300/12)*(0.05/10)*exp(-time(1:64:T)/10.0),'rs','MarkerFaceColor','y')
hold off
