clear
clc


NoOfPoints = 1000;
a  = normrnd(25,30,[NoOfPoints,1]);
b  = normrnd(12,1,[NoOfPoints,1]);

%Condition F1 as non-negative
%Condition F2 to be less than F1

F1 = abs(a);
F2 = abs(b);

F3 = F1+ F2 ;
F4 = F3 -F2 + F1;
%F5 = F4 + F2 - F1;

F = [ F1, F2, F3, F4];

ErrSigma = 0.1;

e1  = normrnd(0,sqrt(ErrSigma),[NoOfPoints,1]);
e2  = normrnd(0,sqrt(ErrSigma*2),[NoOfPoints,1]);
e3  = normrnd(0,sqrt(ErrSigma*3),[NoOfPoints,1]);
e4  = normrnd(0,sqrt(ErrSigma*2),[NoOfPoints,1]);
%e5  = normrnd(0,sqrt(ErrSigma*1),[NoOfPoints,1]);

Z = [ F1+e1, F2+e2, F3+e3, F4+e4];

save InputdataIPCA F Z