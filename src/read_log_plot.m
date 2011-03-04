function handle_data()


a=importdata('1/rlt.dat');
[c r]=size(a)
for i=1:c
    a_left=1.0/a(i,1);
    a_right=1.0/a(i,2);
    c_left=a(i,3);
    c_right=a(i,4);
    b_left=a(i,5);
    b_right=a(i,6);
    x=0:0.001:0.3;
    plot(x,c_left*a_left*(x.^2+a_left^2+b_left^2)./((a_left^2+b_left^2+x.^2).^2-4*x.^2*b_left^2));
    ylim([0,5.0]);
    hold on
    plot(x,c_right*a_right*(x.^2+a_right^2+b_right^2)./((a_right^2+b_right^2+x.^2).^2-4*x.^2*b_right^2));
    hold off
    pause(0.1);
end