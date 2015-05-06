function [ind v] = simplex(A,b,c,m,n,x)
% Devolve -1 se o custo for ilimitado e a direção em que vai para
% infinito é devolvida em v. 
% Devolve 0 se houver solução ótima em ind e em v devolve o a
% solução ótima encontrada.

