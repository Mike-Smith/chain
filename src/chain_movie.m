function chain_movie()

a=importdata('tmp_data_for_propagator_real.dat');

s=size(a);
k=s(1);
s=s(2);
axis tight
set(gca,'nextplot','replacechildren');
% Record the movie
for i=1:k/10
    i;
    v=a(i*10,1:2:s);
    scatter(v,zeros(1,s/2),'filled');
%    view(-60,30);
    axis([-1 13 -1 1]);
    F(i) = getframe;
%    pause(0.01);
end
%movie(F,1,10);
movie2avi(F,'vib.avi') ;