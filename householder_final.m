clc
clear

A = input('Enter Square Symmetric Matrix A:');
[m,n] = size(A);
if(m~=n)
    disp('Matrix is not Square!');
    return
end
if (~isequal(A,A'))
    disp('Matrix is not Symmetric!');
    return
end

A_new = householder_transform_func(A);
disp('Below is the new transformed tridiagonal matrix:');
disp(A_new);

bb = zeros(n,1);
cc = zeros(n - 1,1);
for i = 1:n
    bb(i,1) = A_new(i,i);
end
for i = 1:n-1
    cc(i,1) = A_new(i,i+1);
    if(cc(i,1) == 0)
        disp('Sturm Sequence Method Not Valid!');
        return
    end
end
changes1 = 0;
changes2 = 0;
x1 = 0;
x2 = 0;
choice = input('enter 1 if u want to see sign scheme at each value of sturm sequence else enter 2:');
tol = input('enter absolute tolerance:');
while(changes1 == changes2)
    x1 = input('Enter first guess of the eigen value:');
    changes1 = count(x1,bb,cc,choice);
    if(changes1 == -1)
        sprintf('%d is an eigen value',x1)
        return
    end
    x2 = input('Enter second guess of the eigen value:');
    changes2 = count(x2,bb,cc,choice);
    if(changes2 == -1)
        sprintf('%d is an eigen value',x2)
        return
    end
end

while(abs(x1 - x2)>tol)
    x1new = (x1 + x2)/2;
    changes1new = count(x1new,bb,cc,choice);
    if(changes1new == -1)
        sprintf('%d is an eigen value',x1new)
        return
    end
    if(changes1new == changes1)
        x1 = x1new;
    else
        x2 = x1new;
    end
end
disp('One of the guessed eigen values is:')
disp((x1+x2)/2)
function A_new = householder_transform_func(A)
    A_old = A;
    A_new = A;
    [m,n] = size(A);
    I = eye(n);
    for i = 1:n - 2
    
        w = zeros(n,1); %initialization of normal to plane
        s = sqrt(sum(A_old(i,i+1:n).^2));
        if(s == 0)
            continue;
        end
        w(i+1,1) = sqrt( (  1  +  abs( A_old(i,i+1) )/s  )/2 );
    
        for j = i+2:n
            w(j,1) = ( A_old(i,j) * sign(A_old(i,i+1)) )/(2 * s * w(i+1,1));
        end
    
        P = I - 2 * ( w(1:n,1)* (w(1:n,1))') ;
        A_new = P * A_old * P;
    
        A_old = A_new;
    end
    A_new(abs(A_new(:))<5e-14)=0;%to avoid display of some values as *
end

function changes = count(x,bb,cc,choice)%returns number of changes
    changes = 0;
    [n,~] = size(bb);
    signs = zeros(n+1,1);
    signs(1,1) = 1;
    values = zeros(n,1);
    values(1,1) = x - bb(1,1);
    values(2,1) = (x - bb(2,1))*values(1,1) - cc(1,1)*cc(1,1)*1;
    for i = 3:n
        values(i,1) = (x - bb(i,1))*values(i-1,1) - cc(i-1,1)*cc(i-1,1)*values(i-2,1);
    end
    
   
    for i = 1:n
        if i == 1
            if(values(i,1)>0)
                signs(i,1) = 1;
            elseif(values(i,1)<0)
                signs(i,1) = 0;
                changes = changes + 1;
            else
                signs(i,1) = 1;
            end
        else
            if(values(i,1)>0 && signs(i - 1,1)==1)
                signs(i,1) = 1;
            elseif(values(i,1)>0 && signs(i-1,1)==0)
                signs(i,1) = 1;
                changes = changes + 1;
            elseif(values(i,1)<0 && signs(i-1,1)==1)
                signs(i,1) = 0;
                changes = changes + 1;
            elseif(values(i,1)<0 && signs(i-1,1)==0)
                signs(i,1) = 0;
            elseif(values(i,1) == 0)
                signs(i,1) = signs(i-1,1);
            end
        end
    end
    sign_scheme = strings(n + 1,1);
    sign_scheme(1,1) = '+' ;
    if(choice == 1)
        for i = 1:n
            if(signs(i,1) == 1)
                sign_scheme(i+1,1) = '+';
            else
                sign_scheme(i+1,1) = '-';
            end
        end
    end
    if(choice == 1)
        sprintf('Sign Scheme for %d in the sturm sequence is:',x)
        disp(sign_scheme');
        sprintf('number of changes in sign = %d',changes)
    end
    if(values(n,1)== 0)
        changes = -1;
        sprintf('%d is an eigen value.Enter suitable further guesses by rerunning the program.',x)
        return
    end
end