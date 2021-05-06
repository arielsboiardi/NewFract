function I=colorcell2img(F)
% colorcell2img(F) converte la cell array F che in ogni elemento ha una
% tripletta che rappresenta un colore in un'immagine RBG.

[n,m]=size(F);

for kdx=1:n
    for jdx=1:m
        FR(kdx,jdx)=F{kdx,jdx}(1);
        FG(kdx,jdx)=F{kdx,jdx}(2);
        FB(kdx,jdx)=F{kdx,jdx}(3);
    end
end

I=cat(3,FR,FG,FB);
end
