function [r]=Bairstow_roots(p,tol, Nmax)
% [r]=Bairstow_roots(p,tol, Nmax) determina le radici del polinomio p
% utilizzando il metodo di Bairstow.
%
% Input:    p polinomio di cui vogliamo calcolare le radici
%           tol tolleranza per il test residuale nella fattorizzazione con
%           metodo Bairstow
%           Nmax massimo numeor di iterazioni sia nella divisione che nelle
%           scelte delle approssimazioni inizali
% Output:   r vettore delle radici
%           Niter, numero di scelte della prima approssimazione

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n=numel(p);
deg=n-1; % grado del polinomio
r=[];
Niter=0;

%%% Ci libeariamo di qualche caso semplice a priori
% Non è nulla di particolarmente intelligente, ma direi che vale la pena di
% farlo.
if p(n)==0
    if p(n-1)==0 % Annido i due test per farne solo uno a priori
        r=[0, 0];
        p=p(1:n-2);
        deg=deg-2;
    else
        r=0;
        p=p(1:n-1);
        deg=deg-1;
    end
end

%%% Scelta dell'approssimazione inizale
%%%     ===> Si può fare il modo intelligente?
q0=[1 1 1];

%%% Procedo per deflazione con metodo Bairstow fino ad arrivare,
%%% auspicabilmente, ad una equazione lineare.
while deg>=2 && Niter<Nmax
    [q,p1,~]=Bairstow_factor(p,q0,tol,Nmax);
    
    if isempty(q)
        % Questo succede se si annulla i determinante nella function di
        % fattorizzazione con Bairstow. In tal caso proviamo con un'altra
        % scelta di q0 per massimo Nmax volte.
        
        %        ===> Nuovamente: si può fare il modo intelligente?
        % In teoria anche così va bene: il meotdo converge abbastnza bene
        % quindi non dovrebbero esserci grossi problemi. (è vero????)
        
        q0=randi(300,1,3);
        Niter=Niter+1;
        
        % Il comando continue come in Python fa ripetere i ciclo
        % dall'inizio.
        continue
        
        % Se anche dopo aver provato con altri punti inizali non si trova nulla,
        % ci accontentiamo del risultato ottenuto.
        %        ===> Ha senso???
    end
    
    % Posso quindi risolvere l'equazione quadratica con la formula
    % risolutiva per le equazioni di secondo grado e salvare le radici così
    % ottenute.
    rq=root_quadr(q);
    r=[r,rq];
    
    p=p1;
    %%% Procedo per deflazione e cerco le radici del quoziente
    deg=deg-2;
    % Avendo tolto un fattore quadratico il grado si abbassa di 2
end

if Niter==Nmax
    disp('Non è stato possibile determinare tutte le radici :(');
end

if deg==0
    return
end

%%% Nel caso resti un termine lineare aggiungo un ultimo termine, che anche
%%% in questo caso si calcola banalmente.
if deg==1
    r(end+1)=-p(2)/p(1);
end

end
