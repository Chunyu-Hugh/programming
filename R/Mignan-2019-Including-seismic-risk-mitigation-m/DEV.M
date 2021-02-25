作者：满目山河
链接：https://zhuanlan.zhihu.com/p/112471225
来源：知乎
著作权归作者所有。商业转载请联系作者获得授权，非商业转载请注明出处。

clear all;
close all;
clc;
NP=20;                                %种群数量
D=2;                                  %变量的维数
G=100;                                %最大进化代数
F=0.5;                                %变异算子
CR=0.1;                               %交叉算子
Xs=4;                                 %上限
Xx=-4;                                %下限
%赋初值
x=zeros(D,NP);                        %初始种群
v=zeros(D,NP);                        %变异种群
u=zeros(D,NP);                        %选择种群
x=rand(D,NP)*(Xs-Xx)+Xx;              %赋初值
%计算目标函数
for m=1:NP
    Ob(m)=func2(x(:,m));
end
trace(1)=min(Ob);
%差分进化循环
for gen=1:G
    %变异操作
    %r1,r2,r3,m各不相同
    for m=1:NP
        r1=randi([1,NP],1,1);
        while (r1==m)
            r1=randi([1,NP],1,1);
        end
        r2=randi([1,NP],1,1);
        while (r1==r2)||(r2==m)
            r2=randi([1,NP],1,1);
        end
        r3=randi([1,NP],1,1);
        while (r1==r3)||(r2==r3)||(r3==m)
            r3=randi([1,NP],1,1);
        end
         v(:,m)=x(:,r1)+F*(x(:,r2)-x(:,r3));
    end
    %交叉操作
    r=randi([1,D],1,1);
    for n=1:D
        cr=rand(1);
        if(cr<CR)||(n==r)
            u(n,:)=v(n,:);
        else
            u(n,:)=x(n,:);
        end
    end
    %边界处理，边界吸收
    for n=1:D
        for m=1:NP
            if u(n,m)<Xx
                u(n,m)=Xx;
            end
            if u(n,m)>Xs
                u(n,m)=Xs;
            end
        end
    end
    %选择操作
    for m=1:NP
        Ob1(m)=func2(u(:,m));
    end
    for m=1:NP
        if Ob1(m)<Ob(m)
            x(:,m)=u(:,m);
        end
    end
    for m=1:NP
        Ob(m)=func2(x(:,m));
    end
    trace(gen+1)=min(Ob);
    
end
[SortOb,Index]=sort(Ob);
x=x(:,Index);
X=x(:,1);                                      %最优变量
Y=min(Ob);                                     %最优值
%画图
figure
plot(trace);
xlabel('迭代次数')
ylabel('目标函数值')
title('DE目标函数曲线')
%适应度函数
function value=func2(x)
value=3*cos(x(1)*x(2))+x(1)+x(2);
end