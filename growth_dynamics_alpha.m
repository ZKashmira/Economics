% Shooting algorithm
% Study growth model saddle path for different values of alpha (CD power in production function)

delta=.08;beta=.96;lambda=.5;T=200;sigma=1.01;
A=1;
alpha = [0.2, 0.4, 0.6, 0.8];
cols = ['b--', 'g--', 'y--', 'r--'];
clear h
hold on

for j=1:4
    kss=((1/beta-(1-delta))/A/alpha(j))^(1/(alpha(j)-1));
    k0=lambda*kss;
    ksol(1)=k0;
    
    for t=2:T
    kguess(1)=ksol(t-1);
    kmin=ksol(t-1);kmax=kss;
       while abs(kmax-kmin)>.00000015*kss
       kn=.5*(kmin+kmax);
       kguess(2)=kn;
       stop=0;
       i=2;
       while stop < 1
          i=i+1;
             kguess(i)=A*kguess(i-1)^alpha(j)+(1-delta)*kguess(i-1)-...
             (beta*(A*alpha(j)*kguess(i-1)^(alpha(j)-1)+(1-delta)))^(1/sigma)*...
            (A*kguess(i-2)^alpha(j)+(1-delta)*kguess(i-2)-kguess(i-1));
             if kguess(i)<=kguess(i-1), kmin=kn;stop=1;else,kguess(i)=kguess(i);end
             if kguess(i)>kss, kmax=kn;stop=1;else,kguess(i)=kguess(i);end
      end
      end
       ksol(t)=kguess(2);
    end

    c(1:(T-1))=A*ksol(1:(T-1)).^alpha(j)+(1-delta)*ksol(1:(T-1))-ksol(2:T);c(T)=c(T-1);
    k_lom = A * (k.^alpha(j)) - delta * k;
    k_ts = ksol./kss;
    h(j) = plot(k_ts)
    hold on 
       
end
title('Series of k/k*')
hold off
legend(h, '0.2', '0.4', '0.6', '0.8');

