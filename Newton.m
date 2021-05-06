function [z,Niter,zn,Dzn,Reszn]=Newton(f,df,z0,tol,Nmax)
% Newton calcola la radice dell'equazione f=0 con stima iniziale z0
% Input:    f Function che definisce la fuzione di cui trovare gli zeri
%           df Derivata di f
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

% Commento: Diamo in output parecchie cose prendendo spunto dal codice su
%           Quarteroni, "Matematica numerica", in effetti siccome le
%           calcoliamo gi� qui tanto vale farne qualcosa, alla peggio basta
%           non assegnarle.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Niter=0;

fz0=f(z0);
dfz0=df(z0);

zn=[z0];
Dzn=[NaN];

res=abs(f(z0));
Reszn=[res];

esc=1;
while esc>tol && Niter<Nmax
    % Uso una combinazione del test residuale e del test sullo scarto sperando
    % che questa combinazione sia pi� diagnofistica. Magari ci si pu�
    % accontentare un po' nella scelta di tol
    
    if abs(dfz0)<1e-17
        % Verifica utile, infatti se la derivata fosse proprio nulla il metodo
        % non sarebbe significativo. Ma in ogni caso, anche con la derivata
        % molto piccola, i test di arresto non sono indicativi, quindi �
        % inutle stare a fare i conti.
        
        % Rileviamo che questo non � per forza un problema: se la radice �
        % doppia questo accade abbastanza facilmente se siamo vicini alla
        % radice. Esaminando la successione dei residui, degli scarti e
        % delle radici approssimate possiamo capire se il risultato possa
        % essere corretto o meno.
        
        if res<1e-18
            % Se la derivata � molto piccola ma anche la funzione � molto
            % piccola, significa che siamo vicino ad una radice multipla.
            % Se l'approssimazione � buona considero di avere la radice,
            % quindi faccio in modo di uscire dal ciclo senza andare a fare
            % divisioni pericolose. In effetti questo controllo torna col
            % fatto che la funzione tende a zero pi� rapidamente della sua
            % derivata
            esc=0;
        else
            % Se invece il residuo non � plausibilmente piccolo, significa
            % che la derivata si annulla  e basta, e quindi le iterazioni
            % ci spingono all'infinito. Continuare l'iterazione per� �
            % pericoloso.
            disp('Arresto causa annullamento della derivata')
            z=Inf;
            return
        end
    else
        % Se invece ha senso continuare le iterazioni, costruiamo una
        % approssimazione successiva e aggiorniamo i dati e i parametri d i
        % arresto.
        
        z=z0-fz0/dfz0; % Nuova approssimazione
        
        zn=[zn;z]; % Aggiorno i dati di uscita
        scrt=abs(z-z0);
        Dzn=[Dzn; scrt];
        fz=f(z);
        res=abs(fz);
        Reszn=[Reszn; res];
        
        z0=z; % Riassegno i parametri iniziali
        fz0=fz;
        dfz0=df(z0);
        
        Niter=Niter+1; % Aggiorno i test di arresto
        esc=max(scrt,res); 
    end
end

if Niter==Nmax
    % Se l'algoritmo termina perch� � stato superato il numero massimo di
    % iterazioni senza soddisfare i criteri di stop la rdice potrebbe non
    % essere giusta.
    disp('Raggiunto numero massimo di iterazioni')
end
end
