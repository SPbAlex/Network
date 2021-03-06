clear;
m=1;
n=10000;
t=0:1:(n-1);
t1=0.7;
t2=0.3;
MU=4;
R = exprnd(MU,m,n);
RN=poissrnd(MU,m,n);

k=0;
k1=0;
k2=0;
k3=0;
dx = 0.01; 

figure
hist(R, n)
grid on
X=min(R):(max(R)-min(R))/n:max(R);
disp('max number:');
disp(max(R));
f= exppdf(X, MU);
ff=f*n*((max(R)-min(R))/n);
 hold on
 plot(X,ff,'r')
hold off
figure
hist(RN, n)
 
for i=1:1:n
k=k+R(1,i);
k1=k1+RN(1,i);
end
disp('mat.expect.:');
disp(k/n);
disp('disp:');
disp(var(R));

disp('mat.expect. Puasson:');
disp(k1/n);
disp('disp Puasson:');
disp(var(RN));
 
n1=1;
nr1=1;
k=0;
k1=0;
R1(1,10000)=1;
R1(~R1)=[];
RN1(1,10000)=1;
RN1(~RN1)=[];
for i=1:1:n
if (R(1,i)-t1)>0
    R1(1,n1)=R(1,i)-t1;
    k=k+R1(1,n1);
    n1=n1+1;
end;
if (RN(1,i)-t1)>0
    RN1(1,nr1)=RN(1,i)-t1;
    k1=k1+RN1(1,nr1);
    nr1=nr1+1;
end;
end
disp('mat.expect. step 1:');
disp(k/n1);
disp('number of elements on 2nd iteration:');
disp(n1);
disp('2nd iteration disp:');
disp(var(R1));
disp('mat.expect. step 1 Puasson:');
disp(k1/nr1);
disp('number of elements on second iteration Puasson:');
disp(nr1);
disp('disp on second iteration Puasson:');
disp(var(RN1));

disp('max number of elements');
counts(100)=1;
counts(~counts)=[];
counts=(hist(R1,n1));

h=max(counts);
disp(h);
countOfInterval = ( max(R1) - min(R1) ) / dx; 
x = zeros(length(counts)); 
s = 0; 
for i=1:length(counts) 
counts(i) = counts(i) / n1 / dx; 
x(i) = min(R1) + dx * (i - 1) + dx / 2; 
s = s + counts(i) * dx; 
end 
disp( 's == 1? ' ); 
disp( s ); 
figure
plot( x ,counts );
x(~x)=[];
figure
hist(R1, n1)
grid on
X=min(R1):(max(R1)-min(R1))/n1:max(R1);
f= exppdf(X, MU);
ff=f*n1*((max(R1)-min(R1))/n1);
 hold on
 plot(X,ff,'r')
hold off
figure
hist(RN1, nr1)
 
k=0;
R2(1,10000)=1;
R2(~R2)=[];
n2=1;
for i=1:1:(n1-1)
if (R1(1,i)-t2)>0
    R2(1,n2)=R1(1,i)-t2;
    k=k+R2(1,n2);
    n2=n2+1;
end;
end
disp('mat.expect. step 3:');
disp(k/n2);
disp('number of elements step 3:');
disp(n2);
disp('disp step 3:');
disp(var(R2));

figure
hist(R2, n2)
grid on
X=min(R2):(max(R2)-min(R2))/n2:max(R2);
f= exppdf(X, MU);
ff=f*n2*((max(R2)-min(R2))/n2);
 hold on
 plot(X,ff,'r')
hold off
