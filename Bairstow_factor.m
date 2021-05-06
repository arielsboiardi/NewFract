function [q,p1,Niter]=Bairstow_factor(p,q0,tol,Nmax)
% [q,p1,Niter]=Bairstow_factor(p,q0,tol,Nmax) determina la fattorizzazione
% di p(x) con un fattore qudratico q:
%       p(x) = q(x)*p1(x)
%
% I coefficienti del polinomio in p sono ordinati dal termine di grado
% massimo al termine di grado minimo
%       p(x) = p(n)*x^n + ... + p(2)*x + p(1)
% e i coefficenti in q sono tali che
%       q(x) = x^2 + q(2)*x + q(1)
%
% Input:    p coefficienti del polinomio da fattorizzare
%           q0 prima approssimazione del fattore quadratico
%           tol tolleranza nel test di arresto residuale
%           Nmax massimo numeor di iterazioni
% Output:   q fattore quadragico di p prossimo a q0
%           p1 quoziente della divisione polinomiale
%           Niter numero di iterazioni in cui si è raggiunto il risultato

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Controlli sui dati
% è sempre meglio controllare che i dati siano sensati per il problema che 
% intendiamo risolvere.
if ~isa(p,'double')
    error('p deve essere un vettore numerico di coefficienti')
end
if ~isa(q0,'double')
    error('q0 deve essere un vettore numerico di coefficienti')
end

%%% Casi banali
n=numel(p); % Calcolo a parte per evitare di farlo due volte.
if n<=3
%     disp('Il polinomio è già ridotto ad un solo fattore quadratico')
    q=p;
    p1=p;
    Niter=0;
    return
end
m=numel(q0);
if m>3
    disp('q0 dovrebbe essere quadratico, trascuro i termini di troppo')
    q0=q0(m-2:m);
end

%%% Imposto i dati iniziali
n=n-1; % Grado del polinomio

if q0(1)~=1
    q0=q0/q0(1); % Considero il fattore q0 monico
end

u0=q0(2); % La prima approssimazione del fattore è x^2 + u0*x + v0
v0=q0(3);

%%% Determino il primo quoziente ricorsivamente
% Abbimo considerato la scomposizione
%   p(x) = q(x)*[b(3)*x^(n-2)+...+b(n)*x + b(n+1)] + b(n+2)*(x+p) + b(n+3)
b=zeros(n+3,1);
for kdx=3:n+3
    b(kdx)=p(kdx-2)-v0*b(kdx-2)-u0*b(kdx-1);
end

% Valuto il residuo da questa prima divisione, magari è andata bene al
% primo tentativo

R=p(n)-u0*b(n+1)-v0*b(n);
S=b(n+3)+b(n+2)*u0;

res=norm([R,S],2);

% Se non abbiamo ancora raggiunto la precisione richiesta procediamo con
% una versione semplificata dell'algoritmo di Newton per sistemi

Niter=0;
while res>tol && Niter<Nmax
    %%% Determino il secondo quoziente ricorsivamente
    % Con un seconda divisione di polinomi possiamo valutare le derivate di
    % R ed S come mostrato dalla teoria.
    d=zeros(n+1,1);
    for kdx=3:n+1
        d(kdx)=-b(kdx)-u0*d(kdx-1)-v0*d(kdx-2);
    end
    
    %%% Determino una approssimazione successiva (u,v) con il metodo di Newton
    % Il metodo di Newton per sistemi che viene usato qui non è quello diretto
    % perché abbiamo usato varie derivazioni successive dell'equazione iniziale
    % per determinare le derivate in modo numerico, senza quindi dover poi
    % ricorrere alla differenziazione a questo punto.
    
    D=d(n+1)^2+u0*d(n)*d(n+1)+v0*d(n)^2;
    
    if abs(D)<1e-18
        % In questo caso il sistema non si può risolvere, tanto me con
        % Newton...
        disp('Esco per annullamento del determinante Jacobiano, prova a cambiare q0');
        q=[];
        p1=p;
        return
    end
    
    u=u0-(1/D)*((d(n+1)+u0*d(n))*b(n+2)-d(n)*(b(n+3)+u0*b(n+2)));
    v=v0-(1/D)*(v0*d(n)*b(n+2)+d(n+1)*(b(n+3)+u0*b(n+2)));
    
    u0=u; v0=v; % Aggiorno i dati iniziali
    Niter=Niter+1;
    
    %%% Con questa nuova approssimazione divido il polinomio
    b=zeros(n+3,1);
    for kdx=3:n+3
        b(kdx)=p(kdx-2)-v0*b(kdx-2)-u0*b(kdx-1);
    end

    % Valuto il residuo con questa approssimazione
    R=p(n)-u0*b(n+1)-v0*b(n);
    S=b(n+3)+b(n+2)*u0;
    res=norm([R,S],2);
end

%%% Assegno i due fattori che sono stati determinati
q=[1 u0 v0];
p1=b(3:n+1)';
return