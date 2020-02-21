% Study the saddle path of the growth model for different values of sigma (CRRA)

%Parameter values
delta=.08;alpha=1/3;beta=.96;lambda=.5;T=200;
A=1;
sigma = [0.01, 1.01, 2, 3, 4, 10];
cols = ['b--', 'g--', 'y--', 'r--', 'k--', 'c--'];


k = [0 : 0.01 : 5];
k_lom = A * (k.^alpha) - delta * k;
plot(k, k_lom)
title('Saddle path - Neoclassical Growth Model')
hold on

for j = 1:6
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
             kguess(i)=A*kguess(i-1)^alpha+(1-delta)*kguess(i-1)-...
             (beta*(A*alpha*kguess(i-1)^(alpha-1)+(1-delta)))^(1/sigma(j))*...
            (A*kguess(i-2)^alpha+(1-delta)*kguess(i-2)-kguess(i-1));
             if kguess(i)<=kguess(i-1), kmin=kn;stop=1;else,kguess(i)=kguess(i);end
             if kguess(i)>kss, kmax=kn;stop=1;else,kguess(i)=kguess(i);end
      end
      end
       ksol(t)=kguess(2);
    end

    c(1:(T-1))=A*ksol(1:(T-1)).^alpha+(1-delta)*ksol(1:(T-1))-ksol(2:T);c(T)=c(T-1);
    h(j)=plot(ksol, c, cols(j))       
end

hold off
legend(h, '0.01', '1.01', '2', '3', '4', '10');
