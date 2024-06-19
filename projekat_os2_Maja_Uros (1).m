clear all; close all; clc;

    R = 32; 
    yd = 90; 
    wp = 0.5*pi/R; 
    Kmax = 10; 
    

    Rprime = factor(R); % izracunavanje prostih delilaca broja R
    M = length(Rprime); % izracunavanje broja prostih delilaca broja R

   
    P = unique(perms(Rprime),'rows');  % generisanje liste P koja sadrzi sve razlicite permutacije liste Rprime
    [Np,~] = size(P); % izracunavanje broja razlicitih permutacija liste Rprime

    
    C = nchoosek(0:Kmax, M); % generisanje matrice C koja sadrzi sve kombinacije M-torki iz skupa {0,1,...,Kmax}
    Nc = length(C); % izracunavanje broja kombinacija M-torki iz skupa {0,1,...,Kmax}

    %wq = linspace(0, wp, yd); % kreiranje vektora koji sadrzi yd jednako udaljenih tacaka, pocev od 0 do wp
    br_tac_seg=16;
    deltawtmp=0.5/R/br_tac_seg;
    wtmp=(0:deltawtmp:1)*pi;
    indwq=0;
    for br_1=1:length(wtmp)
        for br_2=1:R/2
            if wtmp(br_1)>2*br_2*pi/R-wp & wtmp(br_1)<2*br_2*pi/R+wp
                indwq=indwq+1;
                wq(indwq)=wtmp(br_1);
            end
        end
    end
    Q = length(wq); 

    deltaopt = 1/yd;
    Aopt = inf;

    % Step 7: for p=1 to Np
    for p = 1:Np
        A = zeros(Nc, 1); % pravljenje vektora A sa Nc redova i jednom kolonom ciji su svi elementi 0
        for c = 1:Nc
            K = C(c, :); % upisivanje c-tog reda iz matrice C u vektor K
            Rp = P(p, :); % upisivanje p-tog reda iz matrice P u vektor Rp
            h = 0;
            I = 0;
            O = 0;
            for m = 1:M-1
                pp = abs(K(m)-K(m+1));
                I = 1;
                for k = (m+1):M
                    I = I*(Rp(k));
                end
                O = O + pp*I;
                
            end
            A(c) = K(1)*R + O + K(M); % izracunavanje broja operacija A po formuli
        end

       
        [A, idx] = sort(A); %sortiranje lista A i C u rastucem redosledu
        C = C(idx, :); 

        for c = 1:Nc
            
                
                if A(c) > Aopt
                    break;
                end
        
                delta = max(abs(freqz(C(c, :), Rp, wq))) %delta=max{|H(exp(jw1))|, |H(exp(jw2))|,..., |H(exp(jwQ))|}

                if delta <= deltaopt
                    deltaopt = delta;
                    Aopt = A(c);
                    K = C(c, :);
                    Rm = Rp;
                end
            
        end
    end

    % Izracunavanje Lk, k=0,1,...M, po formuli
    Lk = zeros(M, 1);
    Lk(1) = -K(1);
    Lk(M) = K(M);
    for k=2:M-1
            Lk(k)= K(k)-K(k+1);
    end
