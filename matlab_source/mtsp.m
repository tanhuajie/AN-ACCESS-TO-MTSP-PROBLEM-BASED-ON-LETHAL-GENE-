tic
clc,clear
%%算法参数
w = 100;
g = 1000;
Pm = 0.1;
Pc = 0.8;
n = 5;
%%经纬度坐标
x = [120.7015202
120.6987175
120.6997954
120.70691
120.7056165
120.7031731
120.6928965
120.6943337
120.6973521
120.6962022
120.7011609
120.6939026
120.6983582
120.7025263
120.6914592
120.6960585
120.7005141
120.6998673
120.6928965
120.6964897
120.6969209
120.7052571
120.7088504
120.7087066
120.7130185
120.6896626
120.6937588
120.6993643
120.7129466
120.7002266
120.7015202
];
y  =[36.37423 
36.37457569
36.37591239
36.37579616
36.37248342
36.37753964
36.37800457
36.37521499
36.37876006
36.37643544
36.37905063
36.38021291
36.38056159
36.38120084
36.38201441
36.38247931
36.38276987
36.37079794
36.37079794
36.36824059
36.37143727
36.36899618
36.37021674
36.36731063
36.36829872
36.36661314
36.36242812
36.38741865
36.37201847
36.36428816
36.37423
];
data = [x';y'];
length = size(data,2);
DATA_PATH = [x',repmat(x(length),1,n);y',repmat(y(length),1,n)];%%用于路径描述
%%距离矩阵
city = zeros(length);
for i = 1:length
    for j = 1:length
        city(i,j) = 6370*acos(cosd(data(1,i)-data(1,j))*cosd(data(2,i))*cosd(data(2,j))+sind(data(2,i))*sind(data(2,j)));
    end
end
%%致死特征
MAX = max(max(city));
death = -(length+n)*MAX;
DEATH = repmat(death,1,n);
%%致死片段分布库
cities = zeros(length+n);
for i = 1:length+n
    if i == 1
       cities(i,:) = [city(1,1:length-1),DEATH,city(1,length)];
    end
    if (1<i)&&(i<length)
        cities(i,:) = [city(i,1:length-1),repmat(city(i,1),1,n+1)];
    end
    if (length<=i)&&(i<length+n)
        cities(i,:) = [death,city(1,2:length-1),DEATH,death];
    end
    if i == length+n
        cities(i,:) = [city(1,1:length-1),DEATH,city(1,length)];
    end
end
%%改良圈算法获得较好父本
Parent = zeros(w,length+n);
for k = 1:w
    %%保证父本没有致死片段
    point = 0;
    while point == 0
        dis = 0;
        ran = randperm(length+n-2);
        route = [1,ran+1,length+n];
        for q = 1:length+n-1
            dis = dis + cities(route(q),route(q+1));
        end
        if dis>0
            point = 1;
        end
    end
    signal = 1;
    while signal == 1
        signal = 0;
        for i = 1:length+n-3
            for j = i+3:length+n
               %%保证无致死片段，且符合改良条件
               a = cities(route(i),route(j-1))+cities(route(i+1),route(j));
               b = cities(route(i),route(i+1))+cities(route(j-1),route(j));
                if (a>0)&&(a<b)
                    route(i+1:j-1) = route(j-1:-1:i+1);
                    signal = 1;
                end
            end
        end       
    end
    Parent(k,route) = 1:length+n;%%编码
end
Parent = Parent/(length+n);
Parent(:,1) = 0;
Parent(:,length+n) = 1;
A = Parent;

%%遗传算法
for k = 1:g
    %%杂交
    B = A;
    c = randperm(w);
    for i = 1:2:w*Pc
        F = 2+floor((length+n-2)*rand(1));
        temp = B(c(i),F:length+n);
        B(c(i),F:length+n) = B(c(i+1),F:length+n);
        B(c(i+1),F:length+n) = temp;
    end
    %%变异
    K = 1+floor(w*rand(1,w*Pm));
    C = A(K,:);
    L = size(K,2);
    for j = 1:L
        bw = 2+floor((length+n-2)*rand(1,3));
        bw = sort(bw);
        C(j,:) = C(j,[1:bw(1)-1,bw(2)+1:bw(3),bw(1):bw(2),bw(3)+1:length+n]);
    end
    %%种群
    G = [A;B;C];
    TL = size(G,1);
    [DX,IX] = sort(G,2);%%解码
    %%淘汰含致死片段的子代
    Fdis = zeros(1,TL);
    for j = 1:TL
        for i = 1:length+n-1
            Fdis(j) = Fdis(j)+cities(IX(j,i),IX(j,i+1));
        end
    end
    pass = find(Fdis>0);
    PASS = IX(pass,:);
    TLL = size(PASS,1);
    %%改良子代
    for m = w+1:TLL
        origin = PASS(m,:);
        signal = 1;
        while signal == 1
            signal = 0;
            for i = 1:length+n-3
                for j = i+3:length+n
                     %%保证无致死片段，且符合改良条件
                     a = cities(origin(i),origin(j-1))+cities(origin(i+1),origin(j));
                     b = cities(origin(i),origin(i+1))+cities(origin(j-1),origin(j));
                    if (a>0)&&(a<b)
                        origin(i+1:j-1) = origin(j-1:-1:i+1);
                        signal = 1;
                    end
                end
            end
        end
        PASS(m,:) = origin;
    end
    %%选取下一次迭代的父本
    show = zeros(1,TLL);
    for j = 1:TLL
        for i = 1:length+n-1
            show(j) = show(j)+cities(PASS(j,i),PASS(j,i+1));
        end
    end
    [DY,IY] = sort(show);
    D = PASS(IY(1:w),:);
    %%父本编码
    H = zeros(w,length+n);
    for i = 1:w
        H(i,D(i,:)) = 1:length+n;
    end
    A = H/(length+n);
    A(:,1) = 0;
    A(:,length+n) = 1;
end
%%迭代后取最优
path = D(1,:);
long = DY(1);
fprintf('最短路径是 %f km    ',long)
toc
TIME = find(path>(length-1));
TLLL = size(TIME,2);
figure(1)
clf
hold on
for i = 1:TLLL
    if i == 1
        newline;
        fprintf('第%d段路径是',i);
        disp(path(1,2:TIME(1)-1)-1);
        newline;
        xx = DATA_PATH(1,path(1,1:TIME(1)));
        yy = DATA_PATH(2,path(1,1:TIME(1)));
        plot(xx,yy,'-o')
    else
        fprintf('第%d段路径是',i);
        disp(path(1,TIME(i-1)+1:TIME(i)-1)-1);
        newline;
        xx = DATA_PATH(1,path(1,TIME(i-1):TIME(i)));
        yy = DATA_PATH(2,path(1,TIME(i-1):TIME(i)));
        plot(xx,yy,'-o')
    end
end
hold off
        
