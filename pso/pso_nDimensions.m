%% ��ջ���
clc;clear all;close all;
%% Ŀ�꺯��
%  ���������������Schaffer��������Сֵ�ڣ�0,0������Ϊ0
fun= @(a,b)(0.5+(sin(sqrt(a.^2+b.^2)).^2-0.5)./((1+0.001.*(a.^2+b.^2)).^2));
%  ��ͼ������������
figure(1);
[x0_1, x0_2]=meshgrid(-5:0.1:5);
y0=fun(x0_1,x0_2);
mesh(x0_1, x0_2, y0);
hold on;
%% ������Ⱥ����
%   ��Ҫ��������
sizepop = 500;                   % ��ʼ��Ⱥ����
dim = 2;                         % �ռ�ά��
ger = 300;                       % ����������     
xlimit = [ -5,5 ; -5 ,5 ];       % ����λ�ò�������(�������ʽ���Զ�ά)
vlimit = [ -1.5,1.5 ; -1.5,1.5]; % �����ٶ�����
c_1 = 0.8;                       % ����Ȩ��
c_2 = 0.5;                       % ����ѧϰ����
c_3 = 0.5;                       % Ⱥ��ѧϰ���� 
%% ���ɳ�ʼ��Ⱥ
%  ����������ɳ�ʼ��Ⱥλ��
%  Ȼ��������ɳ�ʼ��Ⱥ�ٶ�
%  Ȼ���ʼ��������ʷ���λ�ã��Լ�������ʷ�����Ӧ��
%  Ȼ���ʼ��Ⱥ����ʷ���λ�ã��Լ�Ⱥ����ʷ�����Ӧ��
%  ��ͼ����������
 for i = 1:dim
    pop_x(i,:) = xlimit(i, 1) + (xlimit(i, 2) - xlimit(i, 1)) * rand(1, sizepop);%��ʼ��Ⱥ��λ��
end    
pop_v = rand(dim, sizepop);                   % ��ʼ��Ⱥ���ٶ�
gbest = pop_x;                                % ÿ���������ʷ���λ��
fitness_gbest = fun(pop_x(1,:),pop_x(2,:));   % ÿ���������ʷ�����Ӧ��
zbest = pop_x(:,1);                           % ��Ⱥ����ʷ���λ��
fitness_zbest = fitness_gbest(1);             % ��Ⱥ����ʷ�����Ӧ��
for j=1:sizepop
    if fitness_gbest(j) < fitness_zbest       % �������Сֵ����Ϊ<; ��������ֵ����Ϊ>; 
        zbest = pop_x(:,j);
        fitness_zbest=fitness_gbest(j);
    end
end

plot3(gbest(1,:),gbest(2,:),fun(gbest(1,:),gbest(2,:)), 'ro');title('��ʼ״̬ͼ');
hold on;
figure(2);
mesh(x0_1, x0_2, y0);
hold on;
plot3(gbest(1,:),gbest(2,:),fun(gbest(1,:),gbest(2,:)), 'ro');
hold on;

%% ����Ⱥ����
%    �����ٶȲ����ٶȽ��б߽紦��    
%    ����λ�ò���λ�ý��б߽紦��
%    ��������Ӧ����
%    ��������Ⱥ��������λ�õ���Ӧ��
%    ����Ӧ���������ʷ�����Ӧ�����Ƚ�
%    ������ʷ�����Ӧ������Ⱥ��ʷ�����Ӧ�����Ƚ�
%    �ٴ�ѭ�������

iter = 1;                        %��������
times = 1;                       %������ʾ��������
record = zeros(ger, 1);          % ��¼��
while iter <= ger
    
    %    �����ٶȲ����ٶȽ��б߽紦�� 
    pop_v =  c_1 * pop_v  + c_2 * rand *(gbest - pop_x) + c_3 * rand *(repmat(zbest, 1, sizepop) - pop_x);% �ٶȸ���
    for i=1:dim 
        for j=1:sizepop
            if  pop_v(i,j)>vlimit(i,2)
                pop_v(i,j)=vlimit(i,2);
            end
            if  pop_v(i,j) < vlimit(i,1)
                pop_v(i,j)=vlimit(i,1);
            end
        end
    end 
    
    %   ����λ�ò���λ�ý��б߽紦��
    pop_x = pop_x + pop_v;  % λ�ø���
    %   �߽�λ�ô���
    for i=1:dim 
        for j=1:sizepop
            if  pop_x(i,j)>xlimit(i,2)
                pop_x(i,j)=xlimit(i,2);
            end
            if  pop_x(i,j) < xlimit(i,1)
                pop_x(i,j)=xlimit(i,1);
            end
        end
    end
    
    %    ����Ӧ����
    for j=1:sizepop
        if rand > 0.85
            i=ceil(dim*rand);
            pop_x(i,j)=xlimit(i, 1) + (xlimit(i, 2) - xlimit(i, 1)) * rand;
        end
    end
    
    %    ��������Ⱥ��������λ�õ���Ӧ��
    fitness_pop = fun(pop_x(1,:),pop_x(2,:)) ; % ��ǰ���и������Ӧ��
    
    %    ����Ӧ���������ʷ�����Ӧ�����Ƚ�
    for j = 1:sizepop      
        if fitness_pop(j) < fitness_gbest(j)       % �������Сֵ����Ϊ<; ��������ֵ����Ϊ>; 
            gbest(:,j) = pop_x(:,j);               % ���¸�����ʷ���λ��            
            fitness_gbest(j) = fitness_pop(j);     % ���¸�����ʷ�����Ӧ��
        end 
    end
    
    %    ������ʷ�����Ӧ������Ⱥ��ʷ�����Ӧ�����Ƚ�
    for j = 1:sizepop  
        if fitness_gbest(j) < fitness_zbest        % �������Сֵ����Ϊ<; ��������ֵ����Ϊ>; 
            zbest = gbest(:,j);                    % ����Ⱥ����ʷ���λ��  
            fitness_zbest=fitness_gbest(j);        % ����Ⱥ����ʷ�����Ӧ��  
        end
    end
    
    record(iter) = fitness_zbest;%���ֵ��¼
    
    if times >= 10        %��ʾ���� ������ȥ
        cla;
        mesh(x0_1, x0_2, y0);
        plot3(pop_x(1,:),pop_x(2,:),fun(pop_x(1,:),pop_x(2,:)), 'ro');title('״̬λ�ñ仯');
        pause(0.5);
        times=0;
    end
    
    iter = iter+1;
    times=times+1;        %��ʾ���� ������ȥ
end
%% ����������

figure(3);plot(record);title('��������')
figure(4);
mesh(x0_1, x0_2, y0);
hold on;
plot3(pop_x(1,:),pop_x(2,:),fun(pop_x(1,:),pop_x(2,:)), 'ro');title('����״̬ͼ');

disp(['����ֵ��',num2str(fitness_zbest)]);
disp('����ȡֵ��');
disp(zbest);