function find_a_b()
% when calculate the situation with r(t)=c*exp(-at)cos(bt) 
% the peak of \hat{r}(omega) is determined by a, b and c.
% find_a_b() calculate a,b from the position of the peak (ll) and the 
% hieght of the peak (hh), see expcos.mw

% The parameters:
% To make the results for chain_expcos roughly comparable to those of chain (with UO noise), I follow some corresponses.
% In expcos noise case, \gamma(\omega)= 2*c*a*(w^2+a^2+b^2)/((a^2+w^2-2*w*b+b^2)*(a^2+w^2+2*w*b+b^2))
% \gamma(t)=c*exp(-a*t)*cos(b*t)
% In UO noise, \gamma(\omega)=2*e/(1+w^2*tau^2)
% \gamma(t)=e/tau*exp(-t/tau)
% 1. integrations of \gamma(\omega) for UO case and expcos case are roughly the same. which means c=e/tau
% 2. b and a are determined by the position and hieght of the peak.
% 3. The width of \gamma(\omega) for UO case is roughly 1/tau=0.025-0.1fs^-1, So for expcos case, I choose the height of the peak so that the width of \gamma(\omega) to be roughly the same. The results are also displayed in find_a_b.m. It's lucky characteristic frequency of the bond is 0.15fs^-1, so the width of \gamma(\omega) is small enough to distinguish omega and 2omega.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% shift peak position, keep height and width
a_left=1.0/9.988889e+01;
b_left=1.500004e-01;
c_left=2*5e-3;
%a,b,c is defined in expcos.html
ll_left=5e-2;
hh_left=25;
L=num2str(ll_left);
H=num2str(hh_left);
eq1=[L,'=(a^2+b^2)^(1/4)*(2*b-(a^2+b^2)^(1/2))^(1/2)'];
eq2=[H,'=a/(2*b*((a^2+b^2)^(1/2)-b))'];
rlt=solve(eq1,eq2);
s=size(rlt.a);
for j=1:s
    a_left=double(rlt.a(j));
    b_left=double(rlt.b(j));
    if(a_left>0 && b_left>0 && isreal(a_left) && isreal(b_left))
        break
    end
end



a_right=0.0;
b_right=0.0;
c_right=2*5e-3;
%a_right,b_right will be calculated by ll_right and hh_right
ll_right=0.05:0.005:0.25;%position of peak of gamma_r(omega)
hh_right=25;%height of peak of gamma_r(omega)

f=fopen('list.dat','w');
for i=ll_right %change the position of the peak, keep the hieght constant
    L=num2str(i);
    H=num2str(hh_right);
    eq1=[L,'=(a^2+b^2)^(1/4)*(2*b-(a^2+b^2)^(1/2))^(1/2)'];
    eq2=[H,'=a/(2*b*((a^2+b^2)^(1/2)-b))'];
    rlt=solve(eq1,eq2);
    s=size(rlt.a);
    for j=1:s
        a_right=double(rlt.a(j));
        b_right=double(rlt.b(j));
        if(a_right>0 && b_right>0 && isreal(a_right) && isreal(b_right))
            fprintf(f,'./ main -no_output %e %e %e %e %e %e\n',1/a_left,1/a_right,c_left,c_right,b_left,b_right);
            break
        end
    end
     x=-0.1:0.001:0.6;
     plot(x,2*5e-2./(1+x.^2*10^2));%the correlation time is roughly 10fs and 40fs for UO noise.
     hold on
     plot(x,c_left*a_left*(x.^2+a_left^2+b_left^2)./((a_left^2+b_left^2+x.^2).^2-4*x.^2*b_left^2));
     ylim([0,0.5]);
     hold on
     plot(x,c_right*a_right*(x.^2+a_right^2+b_right^2)./((a_right^2+b_right^2+x.^2).^2-4*x.^2*b_right^2));
     hold off
     pause(0.1);
end
fclose(f);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% increase heigth, keep integration and position
% a_left=1.0/9.988889e+01;
% b_left=1.500004e-01;
% c_left=0.05;
% 
% a_right=0.0;
% b_right=0.0;
% c_right=0.05;
% ll_right=0.15;    
% hh_right=50:1:150;
% f=fopen('list.dat','w');
% for i=hh_right
%     L=num2str(ll_right);
%     H=num2str(i);
%     eq1=[L,'=(a^2+b^2)^(1/4)*(2*b-(a^2+b^2)^(1/2))^(1/2)'];
%     eq2=[H,'=a/(2*b*((a^2+b^2)^(1/2)-b))'];
%     rlt=solve(eq1,eq2);
%     s=size(rlt.a);
%     for j=1:s
%         a_right=double(rlt.a(j));
%         b_right=double(rlt.b(j));
%         if(a_right>0 && b_right>0 && isreal(a_right) && isreal(b_right))
%             fprintf(f,'1/ main %e %e %e %e %e %e\n',1/a_left,1/a_right,c_left,c_right,b_left,b_right);
%             break
%         end
%     end
%      x=0:0.001:0.3;
%      plot(x,c_left*a_left*(x.^2+a_left^2+b_left^2)./((a_left^2+b_left^2+x.^2).^2-4*x.^2*b_left^2));
%      ylim([0,5.0]);
%      hold on
%      plot(x,c_right*a_right*(x.^2+a_right^2+b_right^2)./((a_right^2+b_right^2+x.^2).^2-4*x.^2*b_right^2));
%      hold off
%      pause(1);
% end
% fclose(f);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
%% keep position and width, increase height ( and so integration)
% a_left=1.0/9.988889e+01;
% b_left=1.500004e-01;
% c_left=0.025;
% 
% a_right=1.0/9.988889e+01;
% b_right=1.500004e-01;
% c_right=0.05;
% height=0.02:0.002:0.1;
% f=fopen('list.dat','w');
% for i=height
%     c_right=i;
%     fprintf(f,'./ main -no_output %e %e %e %e %e %e\n',1/a_left,1/a_right,c_left,c_right,b_left,b_right);
%     x=0:0.001:0.3;
%     plot(x,c_left*a_left*(x.^2+a_left^2+b_left^2)./((a_left^2+b_left^2+x.^2).^2-4*x.^2*b_left^2));
%     ylim([0,5.0]);
%     hold on
%     plot(x,c_right*a_right*(x.^2+a_right^2+b_right^2)./((a_right^2+b_right^2+x.^2).^2-4*x.^2*b_right^2));
%     hold off
%     pause(0.1);
% end
% fclose(f);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% To observe the 2\omega effect, it's better to increase the peak of \gamma(\omega) 
% %%for both ll_right and ll_left (see in find_a_b.m) from the bond frequency
% %%(0.15) to 3 times (0.45). 
% 
% 
% 
% 
% a_right=0.0;
% b_right=0.0;
% c_right=2*5e-3;
% %a_right,b_right will be calculated by ll_right and hh_right
% ll_right=0.10:0.01:0.60;%position of peak of gamma_r(omega)
% hh_right=50;%height of peak of gamma_r(omega)
% 
% f=fopen('list.dat','w');
% for i=ll_right %change the position of the peak, keep the hieght constant
%     L=num2str(i);
%     H=num2str(hh_right);
%     eq1=[L,'=(a^2+b^2)^(1/4)*(2*b-(a^2+b^2)^(1/2))^(1/2)'];
%     eq2=[H,'=a/(2*b*((a^2+b^2)^(1/2)-b))'];
%     rlt=solve(eq1,eq2);
%     s=size(rlt.a);
%     for j=1:s
%         a_right=double(rlt.a(j));
%         b_right=double(rlt.b(j));
%         a_left=a_right;
%         b_left=b_right;        
%         c_left=c_right;
%         if(a_right>0 && b_right>0 && isreal(a_right) && isreal(b_right))
%             fprintf(f,'./yzhou -digested %e %e %e %e %e %e\n',1/a_left,1/a_right,c_left,c_right,b_left,b_right);
%             break
%         end
%     end
%      x=-0.1:0.001:0.7;
%      plot(x,2*5e-2./(1+x.^2*10^2));%the correlation time is roughly 10fs and 40fs for UO noise.
%      hold on
%      plot(x,c_left*a_left*(x.^2+a_left^2+b_left^2)./((a_left^2+b_left^2+x.^2).^2-4*x.^2*b_left^2),'--rs','color','red');
%      ylim([0,0.5]);
%      hold on
%      plot(x,c_right*a_right*(x.^2+a_right^2+b_right^2)./((a_right^2+b_right^2+x.^2).^2-4*x.^2*b_right^2),'color','blue');
%      hold off
%      pause(0.05);
% end
% fclose(f);