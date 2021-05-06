%% Questo prgramma disegna l'insieme di Fatou 
% Il seguente codice è piuttosto adattbile, anche se mancano alcune
% migliorie. Con qualche modifica permette, credo, di tracciare i bacini
% attrazione in campo complesso per qualunque mappa iterativa per il
% calcolo di radici di euqazioni non lineari

close all
clear all
            %%%%% Parametri %%%%%%
Nmax=20; % Massimo numero di iterazioni
tol=1e-9; % Tolleranza per i criteri di arresto. Non è l'errore nella soluzione

            %%%%%% Regione %%%%%%%
% Rappresento la regione di piano complesso [Rez1,Rez2]x[Imz1,Imz2]
Rez1=-2;
Rez2=2;
Imz1=-2;
Imz2=2;

        %%%%% Risoluzione dell'immagine %%%%
% La cosa è piuttosto dispendiosa, quindi impostiamo la risoluzione
% dell'immagine. Per come costruito il programma, il risultato sarà
% un'immagine di pts^2 pixel.
pts=300; 

% Costruisce il piano complesso con la risoluzioen richiesta
x=linspace(Rez1,Rez2,pts);
y=linspace(Imz1,Imz2,pts);
[X,Y]=meshgrid(x,y);
C=X+1i*Y; 

        %%%%%%% Funzione da studiare %%%%%%%%%%%
syms z
f=(z-0.1).*(z+0.1).*(z.^3-1);
molt=2;

        %%%%%%% Da qui fa da solo %%%%%%%%%%%%%%%%%5
df=diff(f);
ddf=diff(df);
f=matlabFunction(f);
df=matlabFunction(df);
ddf=matlabFunction(ddf);

% f= z.^2.*(z.^3-1);
% f= (z-0.1).*z.*(z.^3-1);

% f=(z-0.01).*z.*(z.^3-1).*log(z);

% f=@(z) z.*(z.^3 - 1).*(z - 1/10);
% df=@(z) z.*(z.^3 - 1) + 3*z.^3.*(z - 1/10) + (z.^3 - 1).*(z - 1/10);
% ddf=@(z) 12*z.^2.*(z - 1/10) + 8*z.^3 - 2;


% f=@(z) (z.^2+1);
% df=@(z) 2*z;

% f=@(z) z.^3-2*z+2;
% df=@(z) 3*z.^2-2;

% f=@(z) 35*z.^9 - 180*z.^7 + 378*z.^5 - 420*z.^3 + 315*z;
% df=@(z) 315*z.^8 - 1260*z.^6 + 1890*z.^4 - 1260*z.^2 + 315;
% ddf=@(z) 2520*z.^7 - 7560*z.^5 + 7560*z.^3 - 2520*z;


Radici=[Inf];
Colori={
    [0 0 0];
    [1 1 1];
    [0.2 0.2 0.8];
    [0.2 0.8 0.2];
    [0.8 0.2 0.2]
    [0.5 0.5 0.7];
    [0.7  0.5 0.5];
    [0.5 0.7 0.7];
    [0.6 0.2 0.8];
    [0.2 0.8 0]};

for idRe=1:pts
    for idIm=1:pts
        z0=C(idRe,idIm); % Punto iniziale
        
        %%% ==> Applico uno dei tre algoritmi <== %%%
        [z,Niter]=Newton_multiroot(f,df,ddf,z0,tol,Nmax);
%         [z,Niter]=Newton(f,df,z0,tol,Nmax);
%         [z,Niter]=Newton_relaxed(f,df,molt,z0,tol,Nmax);

        % Questo programma si può usare con qualunque algoritmo iterativo
        % per il calcolo di radici (più o meno)
        
        if Niter>=Nmax 
            % Scarto le radici che non sono arrivate a convergenza nel
            % numero di iterazioni massimo consentito.
            
            color_z0=[0 0 0];
            F{idRe,idIm}=color_z0;
        else
            % Se la radice è buona guardo quale è la più vicina delle radici
            % già note e cerco di capire se quell appena trovata è nuova oppure
            % cocide ragionevolmente con una di quelle già trovate.
            
            [Dist_radice, idRadice_prox]=min(abs(z*ones(size(Radici))-Radici));
            color_prox_z0=Colori{idRadice_prox};
            
            if Dist_radice<1e-5
                % Se la z trovata con l'algoritmo di Newton a partire da z0
                % è molto vicina ad una radice r già nota, marchiamo z0
                % come facente parte del bacino di attrazione di r.
                
                                %%%%%%%%% Problema (RISOLTO?) %%%%%%%%
                % Problematico capire come definire qui il bacino di attrazione
                % a partire dalla tolleranza sui criteri di stop.
                % Chiaramento vogliamo avere il minor numero possibile di
                % falsi negativi: richieremmo infatti di trovare più bacini
                % delle effettive radici dell'equazione. D'altra parte se
                % due radici sono vicine non vogliamo confonderle!
                
                color_z0=color_prox_z0; % Determino il colore 
                
                % Voglio poi che il colore assegnato a z0 rifletta anche il
                % numero di iterazioni richiesto per raggiungere
                % la'pprossimazione richiesta. Aggiungo quindi al colore
                % assegnato un velo di grigio proporzionale a Niter/Nmax
                color_z0=color_z0+grigio_relativo(color_z0,Niter,Nmax);
                
                F{idRe,idIm}=color_z0; % Salvo il colore nell'immagine
            else
                % Se la radice è nuova, la aggiungo al catalogo delle
                % radici e scelgo un nuovo colore per il suo bacino di
                % attrazione. 
                
                Radici=[Radici;z];
                
                NC=Colori{idRadice_prox+1};
                
                        %%%%%%%% Da fare %%%%%%
                % Sarebbe carino creare una function che scelga i colori da
                % sola. Poi si aggiunge il nuovo colore in fondo all'elenco
                % dei colori in modo da poterlo usare. 
                
                % NC=compl_color(color_prox_z0);
                % Colori=[Colori; NC];
                
                F{idRe,idIm}=NC; % Salvo il colore nell'immagine
            end
        end
    end
end

I=colorcell2img(F);
imshow(I)

