clear all;
% global H
% global alpha;
% global msef;
N=10;
L=20;

% if isempty(H)
%     s=[];
%     s.nel   = 1;
%     s.idx   = int8(zeros(N*L,2));
%     s.h     = complex(zeros(N*L,1),zeros(N*L,1));
%     H       = repmat(s,N,1);
%     alpha=repmat(0.001,N,1); % Set the updating constant
% end
msef=3000;
x=complex(zeros(N,L),zeros(N,L));
d=complex(zeros(N,1),zeros(N,1));
% scope=timescope('AxesScaling','auto','TimeSpanSource','property','TimeSpan',1000);

%SetFull(1,3);
N1=3;
L1=18;
D=round(L1/2);
type = "Full";
obj = AdaptiveMMSE2DS(N,L,N1,L1,2,type);
wght = complex(zeros(N*L,1),zeros(N*L,1));

% SetCross(N1,L1,2); % Set a cross set of coefficients
%SetFull(N1,L1); % Set a cross set of coefficients

del=dsp.Delay(L1/2);
for n=1:N
    fir{n} = dsp.FIRFilter(exp(-(0:1:3)));
end

for SNR=(-6:3:30)
    clear  AdaptMMSE2D
    sigma=10^(-SNR/20.)/sqrt(2);
    for k=0:3000
        
        p=mod(k,L)+1;
        x(:,p) = complex(randn(N,1),randn(N,1))/sqrt(2);
        
        
        % Reference
        d=del(x(:,p)');

      
        % Channel and noise
        for n=1:N
            x(n,p)=fir{n}(x(n,p));
        end
        x(:,p) = x(:,p) + sigma*complex(randn(N,1),randn(N,1));

        % MMSE
        
        
        [y,mse,pow,ww]= obj.Equalizer(x,p,d',wght);
     %   [y,mse,pow]= AdaptMMSE2D(x,p,d');
        %fprintf("%d\t%f\t%f\t%e\n",k,x(1,1+mod(p-1-D,L)),y(1),mse);
        %scope(-10*log10(mse),10*log10(pow));
        wght = ww;
        
    end
  
    fprintf("%d\t%5.1f\t%5.1f\t%5.1f\n",k*N,SNR,-10*log10(mse),10*log10(pow));
    figure(1);
    for i = 1:55
    oo(i,1) = wght(i,1);
    end
    
   % subplot(2,1,1);
    plot(real(oo),'-');
    xticks((0:1:55))
    grid on;
    ylim([-1 1])
    hold on;
    
%     subplot(2,1,2);
%     plot((0:H(2).nel-1),real(H(2).h(1:H(2).nel)/pow));
%     xticks((0:1:30))
%     ylim([-1 1]);
%     grid on;
%     hold on;
end

