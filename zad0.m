%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      maksimizacija funkcije                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
clc

%% Inicijalizacija

% definisanje parametara GA
N=50; L=24; pc=0.9; pm=0.001; G=0.8;%L format prvih 12 bita su x ostali su y 

% crtanje funkcije
fun = @(x,y) (x+1).^2.*(cos(pi*x)+sin(x)).*exp(-abs(x))./(1+y.^2);
fsurf(fun,[-10 10 -5 5]); xlabel('x'); ylabel('y')

% inicijalna populacija
for i = 1:N
    tmp = rand;%odredjuje gen
    gen(i,1:L) = dec2bin(round(tmp*(2^L-1)),L);
end
% gen je matrica N x L i cine je biti cele generacije
%%

uslov = 0;
br_gen = 0;
Fmax = [];

while uslov==0
    
    br_gen = br_gen+1;
       
    % evaluacija generacije
    x = bin2dec(gen(:,1:L/2))/(2^(L/2)-1)*20-10;  %konverzija gen u x u tekstu je x...
    % izmedju -10 i 10 pa treba preskaliratibrojeve koji su od 0 do 1
    y = bin2dec(gen(:,L/2+1:L))/(2^(L/2)-1)*10-5; % konverzija gen u y
    % izmedju -5 i 5 pa treba preskaliratibrojeve koji su od 0 do 1
    
    %f = (x+1).^2.*(cos(pi*x)+sin(x)).*exp(-abs(x))./(1+y.^2)+1;% izracuna
    f= fun(x,y)+1;%jednostavnije i brze +1 je samo da bi uvek bilo pozitivno trayimo max
    %f=-fun(x,y)+3;% Kad bi trazili minimum
    
    [fmax,imax] = max(f);
    Fmax = [Fmax fmax];% Pamti napredak max
    
    % prikaz rezultata
    close all
    
    figure
    fsurf(fun,[-10 10 -5 5]); xlabel('x'); ylabel('y')
    hold all
    plot3(x,y,fun(x,y),'k*','LineWidth',2)
    plot3(x(imax),y(imax),fun(x(imax),y(imax)),'r*','LineWidth',2)
    title((['Generacija broj ' num2str(br_gen) ', najbolji rezultat: ' num2str(max(f))]));
    hold off
    view([0 90]);
    
    pause
    
    % formiranje nove generacije
    gen_staro = gen;
    clear gen
    
    N_reprodukcija = round((1-G)*N);
    N_ukrstanje = N-round((1-G)*N);
    if mod(N_ukrstanje,2)==1
        N_ukrstanje = N_ukrstanje+1;
        N_reprodukcija = N_reprodukcija-1;
    end
    
    % reprodukcija tockom ruleta
    cc = cumsum(f); %kumulativna suma
    [max_v,max_i]=max(f);
    gen(1,1:L)=gen_staro(max_i,1:L);%prvi hromozom uzimamo najbolji hromozom 
    %od predhodne gen i zelim da ga sacuvam
    
    %vrtimo rulet
    for i = 2:N_reprodukcija
        pom = rand*cc(N); %mesto gde se rulet zaustavi
        pom_i = find(sign(cc-pom)==1,1,'first');% nadji prvu kojoj je cc-pom pozitivna 
        gen(i,1:L) = gen_staro(pom_i,1:L); % uzmi taj gen i prebaci ga u sl gen 
    end
    
    % formiranje novih jedinki ukrstanjem
    for i = 1:N_ukrstanje/2
        
        % izbor roditelja
        pom = rand*cc(N);
        pom_i = find(sign(cc-pom)==1,1,'first');
        roditelj1 = gen_staro(pom_i,1:L);
        
        pom = rand*cc(N);
        pom_i = find(sign(cc-pom)==1,1,'first');
        roditelj2 = gen_staro(pom_i,1:L);
        
        % ukrstanje
        if rand<pc
            tacka_ukrstanja = ceil(rand*(L-1));% random delimo gen i ukrstamo
            gen(N_reprodukcija+2*i-1,1:L) = [roditelj1(1:tacka_ukrstanja) ...
                roditelj2(tacka_ukrstanja+1:L)];
            gen(N_reprodukcija+2*i,1:L) = [roditelj2(1:tacka_ukrstanja) ...
                roditelj1(tacka_ukrstanja+1:L)];
        else
            %ako se ne ukrse samo predju u next gen
            gen(N_reprodukcija+2*i-1,1:L) = roditelj1;
            gen(N_reprodukcija+2*i,1:L) = roditelj2;
        end
        
    end
    
    % mutacija
    %round(pm*N*L) treba da je veci od 1
    for i = 1:round(pm*N*L)
        N_i = ceil(rand*(N));%random gen
        L_i = ceil(rand*(L));%random bit
        gen(N_i,L_i) = num2str(1-str2num(gen(N_i,L_i))); %mutiraj
    end
    
    % uslovi za prestanak algoritma
    fs=sort(f);
    if abs(fs(N)-fs(N-10))<0.0001 || br_gen>50
        uslov=1;
    end
    
end

figure
plot(Fmax)