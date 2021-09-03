classdef AdaptiveMMSE2DS
    
    properties
        input;
        N;
        L;
        N1;
        L1;
        D;
        weight;
        pointer;
        msef=3000;
        Pow=1;
        MSE;
        nav=1000;
        s;
        H;
        alpha;
        activeIdx;
        index;
        type;
        indices;
        CC;
    end
    
    methods
        function obj = AdaptiveMMSE2DS(N,L,N1,L1,D,type)%%constructore
            obj.N = N;
            obj.L = L;
            obj.N1 = N1;
            obj.L1 = L1;
            obj.D = D;
            obj.s = struct('nel', 1, 'idx', int8(zeros(N*L,2)), ...
                'h', complex(zeros(N*L,1),zeros(N*L,1)));
            obj.H = repmat(obj.s,N,1);
            obj.alpha=repmat(0.001,N,1);
            obj.MSE = zeros(N,1);
            obj.Pow = ones(N,1);
            obj.activeIdx = zeros (L,2);
            obj.CC = 1;
           
            
            if type == "Cross"
                
                N2 = floor(N1/2);
                if (D<0 || D>L1)
                    D = 0;
                end
                for n =1:size(obj.H)
                    l=1;
                    for k1=0:L1-1        %time index depend on memory
                        obj.H(n).idx(l,1)=k1;
                        obj.H(n).idx(l,2)= 0;      %for each subcarrier(N) that we choose by n we put index of time independ on memory and index of frequency equal to 0
                        l=l+1;
                    end
                    for n1=0:N1-1        %frequancy index depend on window which is here is 3
                        if(n1==N2)      %the index which is common and choose once in a time domain
                            continue;
                        end
                        obj.H(n).idx(l,1)=D;
                        obj.H(n).idx(l,2)= n1-N2;  %at the time index k0 and get a index of frequency up and down
                        l=l+1;
                    end
                    obj.H(n).nel=l-1;
                end
                
            elseif type == "Full"
                obj.CC = 0;
                N2=floor((obj.N1)/2);
                for n=1:size(obj.H)
                    l=1;
                    for n1=0:(obj.N1)-1
                        for k1=0:(obj.L1)-1
                            obj.H(n).idx(l,1)=k1;
                            obj.H(n).idx(l,2)= n1-N2;
                            l=l+1;
                        end
                    end
                   obj.H(n).nel=l-1;
                end
                
            elseif type == "Rhombus"
                obj.CC = 0;
                N2=round(N1/2);
                L2=round(L1/2);
                for n=1:size(obj.H)
                    l=1;
                    for n1=-N2:N2
                        for k1=-L2:L2
                            %if sqrt((n1/N2)^2+(k1/L2)^2)>1
                            %if max(abs(n1/N2),abs(k1/L2))>1
                            if abs(n1/N2) + abs(k1/L2)> 1
                                continue
                            end
                           obj.H(n).idx(l,1)=k1+L2;  
                           obj.H(n).idx(l,2)= n1;
                            l=l+1;
                        end
                    end
                    obj.H(n).nel=l-1;
                end
                
            end
        end
        
        
        %%%%%taking the indices from "FindIndices" class
        
        
        
        function [output,MSEv,Powv,wieghts]=Equalizer (obj,input,ptr,pilot,wght)
            wght1 = reshape (wght,[],1);
           
            wieghts = zeros(200,10);
            for n = 1:obj.N
                output = complex(0);
                for j = 1:obj.H(n).nel
                    obj.activeIdx(j,1) = mod(ptr-1-obj.H(n).idx(j,1),obj.L)+1;
                    obj.activeIdx(j,2) = mod(n-1-obj.H(n).idx(j,2),obj.N)+1;
                    output = output + input(obj.activeIdx(j,2),obj.activeIdx(j,1)) * wght1(j+((n-1)*obj.H(n).nel*obj.CC));
                end
                y(n)=output/obj.Pow(n);  % Useful signal with unitary power
                obj.Pow(n) = obj.Pow(n) + (abs(output)^2 - obj.Pow(n))/obj.nav;
                obj.MSE(n) = obj.MSE(n) + (abs(y(n)-pilot(n))^2 -obj.MSE(n))/obj.nav;
                if(obj.nav < obj.msef)
                    obj.nav=obj.nav+1;
                end
                updateState = obj.alpha(n) *(output - pilot(n));
                
                for j=1:obj.H(n).nel
%                     obj.activeIdx(j,1) = mod(ptr-1-obj.H(n).idx(j,1),obj.L)+1;
%                     obj.activeIdx(j,2) = mod(n-1-obj.H(n).idx(j,2),obj.N)+1;
                    obj.H(n).h(j)= wght1(j+((n-1)*obj.H(n).nel*obj.CC)) - updateState*conj(input(obj.activeIdx(j,2),obj.activeIdx(j,1)));
                end
               wieghts = obj.H(n).h;
            end
            MSEv=obj.MSE;
            Powv=obj.Pow;
        end
    end
    
end


