function s = SisEnacb(x,f,df)
%SISENACB Summary of this function goes here
%   Detailed explanation goes here

%% Init
nx = length(x);
neq = nx-2;

A = zeros(neq);
b = zeros(neq,1);

dfa = df(x(1));
dfb = df(x(end));

fx = f(x);

%% Sistem

% A
for i = 1:neq
    % Prvi clen
    if i ~= 1
        A(i,i-1) = 1/(x(i+1)-x(i+1-1));
    end
    
    % Drugi clen
    A(i,i) = 2 * (1/(x(i+1)-x(i+1-1)) + 1/(x(i+1+1)-x(i+1)));
    
    % Tretji clen
    if i ~= neq
        A(i,i+1) = 1/(x(i+1+1)-x(i+1));
    end
end

% b
for i = 1:neq
    d1 = (fx(i+1)-fx(i+1-1)) / (x(i+1)-x(i+1-1));
    d2 = (fx(i+1+1)-fx(i+1)) / (x(i+1+1)-x(i+1));
    b(i,1) = 3 * (d1/(x(i+1)-x(i+1-1)) + d2/(x(i+1+1)-x(i+1)));
end

b(1,1) = b(1,1) - dfa/(x(2)-x(1));
b(end,1) = b(end,1) - dfb/(x(end)-x(end-1));

%% Solve
s = A\b;
s = [dfa, s', dfb];

end

