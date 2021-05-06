function [z,Niter,zn,Dzn,Reszn]=Newton_multiroot(f,df,ddf,z0,tol,Nmax)
% Newton_multiroot trova radice di f=0 com metodo Newton radici multiple
%
% Un metodo per trovare le radici multiple dell'equazione f=0 senza perdere
% la convergenza quadratica è considerare una funzione f/df che ha le
% stesse radici di f ma le ha smeplici, come provato dalla teoria. In
% questo caso abbiamo bisogno della derivata seconda
%
% Input:    f Function che definisce la fuzione di cui trovare gli zeri
%           df Derivata di f
%           ddf Derivata seconda di f
%           z0 Approssimazione iniziale della radice
%           tol tolleranza per le stime dell'errore (non coincide
%           necessariamente con la tolleranza sull'errore
%           nell'approssimazione finale
%           Nmax Massimo numero di iterazioni consentito
% Output:   z Approssiamzione della radice
%           Niter Numero di iterazioni impegato per raggiungere
%           la'approssimazione restituita
%           zn Successione delle approssimazioni
%           Dzn Differenze fra approssimazioni successive
%           Reszn Residui ad ogni approssimazione
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Niter=0;

fz0=f(z0);
dfz0=df(z0);
ddfz0=ddf(z0);

zn=[z0];
Dzn=[NaN];

res=abs(f(z0));
Reszn=[res];

esc=1;
while esc>tol && Niter<Nmax
    % Uso una combinazione del test residuale e del test sullo scarto sperando
    % che questa combinazione sia più diagnofistica. Magari ci si può
    % accontentare un po' nella scelta di tol
    
    D=(dfz0^2 - fz0*ddfz0);
    N=fz0*dfz0;
    
    if abs(D)<1e-17
        % Se questo denominatore è sotto la precisione macchina la
        % divisione non è significativa. In tal caso però dovremmo avere
        % anche f molto piccola, almeno nel caso di radici doppie, in caso
        % contrario l'aloritmo non converge.
        
        if abs(N)<1e-18
            % Se f è molto piccola possimao considerare di essere finiti su
            % una radice doppia, che ha reso molto piccolo il denominatore
            % senza però aver soddisfatto il criterio dello scarto. Questo
            % può accadere in caso di convergenza molto veloce.
            esc=0;
        else
            % Se invece il residuo non è plausibilmente piccolo, siamo in
            % una situazione di non covnergenza. 
            disp('Arresto causa annullamento del denominatore')
            z=Inf;
            return
        end
    else
        % Se invece ha senso continuare le iterazioni, costruiamo una
        % approssimazione successiva e aggiorniamo i dati e i parametri di
        % arresto.
        
        z=z0-N/D; % Nuova approssimazione
        
        zn=[zn;z]; % Aggiorno i dati di uscita
        scrt=abs(z-z0);
        Dzn=[Dzn; scrt];
        fz=f(z);
        res=abs(fz);
        Reszn=[Reszn; res];
        
        z0=z; % Riassegno i parametri iniziali
        fz0=fz;
        dfz0=df(z0);
        ddfz0=ddf(z0);
        
        Niter=Niter+1; % Aggiorno i test di arresto
        esc=max(scrt,res); 
    end
end

if Niter==Nmax
    % Se l'algoritmo termina perché è stato superato il numero massimo di
    % iterazioni senza soddisfare i criteri di stop la rdice potrebbe non
    % essere giusta.
    disp('Raggiunto numero massimo di iterazioni')
end
end