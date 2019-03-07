%{
Himabshu Rajoria
16IM30008

Cutting stock problem
roll length = 60

required widths and numbers->
55:200
45:150
25:100
35:50
40:175

%}

% An initial basis of pattern that satisfy
B = [ [1 0 0 0 0]
      [0 1 0 0 0]
      [0 0 2 0 0]
      [0 0 0 1 0]
      [0 0 0 0 1] ]

% Calculating initial B_inverse
Binv = inv(B)

% the b array ( required number of rolls)
b = [200;
     150;
     100;
      50;
     175]

% Since any feasible pattern will always cost 1 roll
cost = [1 1 1 1 1]      


%----Delayed column generation method----%

%Objective function
f = -cost * Binv

% Feasibility of pattern
Alp = [ 55 45 25 35 40]
blp = [60]


Aeq =[]
beq = []
lb = [0;0;0;0;0]
ub = []




%since total cuts of a length must be integer
intcons = 1:5

%x is the generated column
%fval is the objective 
[x,fval] = intlinprog(f,intcons,Alp,blp,Aeq,beq,lb,ub)

fval = -fval


% optimality test (zj-cj)
while fval > 1
    
    %--Minimum ratio test--%
    
    bpj = Binv*x
    rhs = Binv*b
    
    % doing the minimum ratio test
    mrt_column = rhs./bpj
    
    % m is theta value
    [m,leaving_column] = min(mrt_column)
    
    % changing the basis
    for row = 1:5
        B(row,leaving_column) = x(row)
    end
    
    Binv = inv(B)
    
    % new zj
    f = -cost * Binv
    
    % finding best possible column
    [x,fval] = intlinprog(f,intcons,Alp,blp,Aeq,beq,lb,ub)
    fval = -fval
end

num_patterns = Binv *b
num_rolls = [1 1 1 1 1]*num_patterns

fprintf('The types of patterns are\n')
B

fprintf('\nthe types of each pattern are \n')
num_patterns

fprintf('The total rolls are ')
num_rolls




    
    
    
 
    
     
        






