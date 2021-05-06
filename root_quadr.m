function rq=root_quadr(q)
% rq=root_quadr(q) calcola le due radici del polinomio quadratico con
% coefficienti q:
%       q(x)=q(1)*x^2+q(2)*x+q(3)

n=numel(q);
if n>3
    warning('q deve avere grado 2, considero solo la parte qiuadratica');
    q=q(n-2:n);
end

D=q(2)^2-4*q(1)*q(3);

rq(1)=(-q(2)+sqrt(D))/2;
rq(2)=(-q(2)-sqrt(D))/2;

end