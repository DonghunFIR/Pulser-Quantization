%% field II _2D linear BeamField _Probe Specification =====================2022/12/13 _WooJin

clc; clear; close all;


addpath('Z:\LAB Common\Z Lab common resources\Software tools\fieldii_20220216\Field_II_ver_3_30_windows')
disp('Beamfield Simulation- Field II')

field_init(0);

%% User define
R                   = 15;
User.txFoci_z       = 15;                 % interest tx focal depth [mm]
longitude           = deg2rad(0);                  % XY plane angle (theta) -> X
latitude            = deg2rad(0);                  % Z  plane angle (pi)       
focusing_range      = 5e-3;                        % Therapeutic range
%% Range angle

xy_range            = 15; 
z_range             = 15;
Data_num            = xy_range * z_range;

%% parameter setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Acoustic Spec. ========================= acous
acous.f0     = 2e6;                          % Center frequency of center frequency [Hz]
acous.FBW    = 1.0;
acous.fs     = 100e6;                        % Sampling frequency
acous.c      = 1540;
acous.lambda = acous.c/acous.f0;             % Wave length [m]
acous.k      = 2*pi/acous.lambda;            % Wave number k  : spatial frequency of a wave
acous.Ts     = 1/acous.fs;                   % Sampling period [s]
acous.acousZ = 1.480e6;                      % Characteristic acoustic impedance [kg/(m^2 s)]
acous.Tprf   = 1/5e3;                        % Pulse repetition frequency [s]


%   Transducer Geometry ==================== trans
% 2D : Symmetric structure
trans.noElement_x   = 8;                    % number of physical elements
trans.noElement_y   = 8;                    % number of physical elements

trans.pitch_x   = 1*1e-3;
trans.pitch_y   = trans.pitch_x;
trans.kerf_x    = 80*1e-6; 
trans.kerf_y    = trans.kerf_x;

trans.width     = trans.pitch_x- trans.kerf_x;
trans.height    = trans.pitch_y - trans.kerf_y;

trans.enabled   = ones(trans.noElement_x, trans.noElement_y);

        
focus_x         = R * sind(longitude);
focus_y         = R * sind(latitude);
focus_z = R * cosd(latitude)*cosd(longitude);

%   Tx Spec. ============================== tx
tx.focus_lens = [focus_x focus_y focus_z]*1e-3;               % Elevational focus [m] : Mechanical focusing
tx.focus_tx   = [focus_x focus_y focus_z]*1e-3;               % transmit focus [m] : Electrical focusing


delay_cell      = zeros(trans.noElement_x *trans.noElement_y  , Data_num);
delay_map_cell  = cell(1,Data_num);
idx             = 1;
delay_map       = zeros(trans.noElement_x,trans.noElement_y);
unique_cell     = cell(1,Data_num);

focus_hisx      = zeros(1,Data_num);
focus_hisy      = zeros(1,Data_num);
focus_hisz      = zeros(1,Data_num);

focus_hisx2     = zeros(z_range,xy_range);
focus_hisy2     = zeros(z_range,xy_range);
focus_hisz2     = zeros(z_range,xy_range);

focus_his       = cell(Data_num , 1);
focus_his_now   = zeros(1,3);
delay_origin    = zeros(trans.noElement_x *trans.noElement_y  , Data_num);

for pi = 1 : z_range
    
    latitude         = pi;                % XY plane angle (theta)
    
    for theta = 1 : xy_range 
               
        longitude        = theta;   
        
        focus_x = R * sind(longitude-1);
        focus_y = R * sind(latitude-1);
        focus_z = R * cosd(latitude-1)*cosd(longitude-1);

        tx.focus_lens = [focus_x focus_y focus_z]*1e-3;
        tx.focus_tx   = [focus_x focus_y focus_z]*1e-3;               % transmit focus [m] : Electrical focusing
        
         
 
        Th2D = xdc_2d_array ( trans.noElement_x, trans.noElement_y, trans.width, trans.height, ...
            trans.kerf_x, trans.kerf_y, trans.enabled, 1, 1, tx.focus_lens );

        delay_each           = xdc_get(Th2D,'focus');
        delay_each           = delay_each(2:end);
        delay_origin( : ,idx)= delay_each;
        delay_each           = round(delay_each,8);
        delay_cell( : , idx) = delay_each;


        for i = 1:trans.noElement_x
            delay_map(i, :) = delay_each( (i-1) * trans.noElement_x + 1 : ...
                i * trans.noElement_x);    
        end

        delay_map_cell{idx} = delay_map;


        focus_hisx(idx) = (focus_x);
        focus_hisy(idx) = (focus_y);
        focus_hisz(idx) = (focus_z);
        focus_his_now(1)   = focus_x;
        focus_his_now(2)   = focus_y;
        focus_his_now(3)   = focus_z;
        focus_his{idx}    = focus_his_now; 
        for i = 1:10
            focus_hisx2(i, :) = focus_hisx( (i-1) * xy_range + 1 : ...
                i * xy_range);    
            focus_hisy2(i, :) = focus_hisy( (i-1) * xy_range + 1 : ...
                i * xy_range); 
            focus_hisz2(i, :) = focus_hisz( (i-1) * xy_range + 1 : ...
                i * xy_range); 
        
        end

        idx = idx + 1;

    end

end

%%% find Numbers of Pulser about each focus
for i = 1 : Data_num
    unique_cell{i} = unique(delay_map_cell{i});
end

unique_cell_stnadard = cell(xy_range,z_range);

idx = 1 ;
col_idx = 1 ;

for i = 1:Data_num

    unique_cell_stnadard{col_idx , idx} = unique_cell{i};
    idx = idx +1;

    if idx == z_range +1 
        idx = 1;
        col_idx = col_idx +1 ;
    end

end

Pulser_num = zeros(xy_range,z_range);
Pulser     = zeros(1,Data_num);

for P = 1:Data_num
    [Pulser(P),~] = size(unique_cell{P});
end

for xy = 1 : xy_range 
    for z = 1 : z_range 
        [Pulser_num(xy,z),~] = size(unique_cell_stnadard{xy,z}); 

    end
end

Pulser_num_trans = zeros(z_range,xy_range);

for i = 1:z_range
    Pulser_num_trans(i,:) = Pulser_num(:,i);
end
%%% find Numbers of Pulser about each focus

% 
% %%% Video 
% figure(2),
% for i = 1:Data_num
% 
%     imagesc(trans.pitch_x*delay_map_cell{i}); grid on; axis image
%     xlim([0.5 8.5])
%     ylim([0.5 8.5])
% 
%     hold on;
%     scatter( (4.5 + round(focus_x)) , 4.5 + round(focus_y),'filled' );
%     
% 
%     hold off;
%     pause(0.00005)
% end
% %%% Video 

 %% PULSER QUANTIZATION   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Pulser_fre       = 80*1e6;                     % Pulser 의 주파수
        CLK_period       = 1 / Pulser_fre;              % CLK 길이
        Practical_Delay     = zeros(64,Data_num);
        times               = zeros(1,1);
        Pulser_level        = zeros(Data_num,1);
        practical_pulser_num= zeros(1,Data_num);
        practical_pulser_num_cell = zeros(xy_range,z_range);
        %%% MAX14808 clk frequency = 160MHz -> 6.25*1e-9

        for i = 1 : Data_num
            Pulser_level(i) = round(max(abs(delay_cell(:,i)))/(CLK_period));
        end

        PM_his = zeros(64,Data_num);

        for u = 1 : Data_num
            for i = 1 : 64
                if delay_cell(i,u) < 0
                    PM_his(i,u) = 1;
                elseif delay_cell(i,u) > 0
                    PM_his(i,u) = 2;
                elseif delay_cell(i,u) == 0
                    PM_his(i,u) = 0 ;
                end
            end
        end

        for n = 1 : Data_num
            now_pulser_level = Pulser_level(n);

            for idx = 1 : now_pulser_level

                for i = 1 : trans.noElement_x * trans.noElement_y

                    if ( 0 < abs(delay_cell(i,n)) ) && ( abs(delay_cell(i,n)) < 0.5 * CLK_period )
                        Practical_Delay(i,n) = 0 ;

                    elseif ( ( 0.5 + idx - 1 ) * CLK_period < abs(delay_cell(i,n)) ) && ( ...
                            abs(delay_cell(i,n)) < ( 0.5 + idx ) * CLK_period )
                        Practical_Delay(i,n) = CLK_period * idx ;

                    end
                end


            end
        end

        for i = 1:64
            if PM_his(i) == 1
                Practical_Delay(i) = Practical_Delay(i) * (-1);
            elseif PM_his(i) == 2
                Practical_Delay(i) = Practical_Delay(i);
            elseif PM_his(i) == 0
                Practical_Delay(i) = Practical_Delay(i);
            end
        end

        Pulser_number = size(unique(Practical_Delay));


%        xdc_focus_times(Th2D,times,Practical_Delay); % Times = zeros(1,1)로 만들어서 넣어야함


        for i = 1:Data_num
            [practical_pulser_num(1,i),~] = size(unique(Practical_Delay(:,i)));
        end
        
        for i = 1:xy_range
            practical_pulser_num_cell(i,:) = practical_pulser_num((i-1) * xy_range + 1 : i * xy_range);    
        end


%% delay 는 8 X 8 이니 64개로 바꿔서 ( ! 밖에서 자르고 넣어야함 ! )
figure(1),
axis image
%%% Pulser numbers
Pulser01                = normalize(Pulser,'range');
delay_origin_diff       = zeros(63,Data_num);
practical_remove        = zeros(1,Data_num);
Pulser_deg              = zeros(z_range , xy_range);
practical_pulser_num01  = practical_pulser_num/(trans.noElement_x*trans.noElement_x);
%practical_pulser_num01  = normalize(practical_pulser_num,'range');
practical_pulser_num01_weight = zeros(1,Data_num);
hold on
grid on
 
% ka=0 ;
% for i = 1 : Data_num
%     if (practical_pulser_num01(i) < 0.80)
%         practical_pulser_num01_weight(i) = practical_pulser_num01(i)/2;
%         ka=ka+1;
%     elseif (practical_pulser_num01(i) > 0.80)
%         practical_pulser_num01_weight(i) = practical_pulser_num01(i);
%     end
% end

practical_pulser_num25 = zeros(1,Data_num);
practical_pulser_num50 = zeros(1,Data_num);
practical_pulser_num75 = zeros(1,Data_num);
practical_pulser_num100 = zeros(1,Data_num);
k50=0;
k75=0;
k100=0;
km=0 ;
for i = 1 : Data_num
    if (practical_pulser_num01(i) < 0.35)
        practical_pulser_num25(i) = practical_pulser_num25(i);
        km=km+1;
    elseif ( practical_pulser_num01(i) >=0.35 && practical_pulser_num01(i) <= 0.50)
        practical_pulser_num50(i) = practical_pulser_num01(i);
        k50 = k50+1;
    elseif ((0.50 <practical_pulser_num01(i)) && (practical_pulser_num01(i) < 0.75))
        practical_pulser_num75(i) = practical_pulser_num01(i);        
        k75 = k75+1;
        
    elseif ((0.75 <=practical_pulser_num01(i)) && (practical_pulser_num01(i) < 1))
        practical_pulser_num100(i) = practical_pulser_num01(i);  
        k100 = k100+1;
        
    end
end
%%
rgb1=[0 1 0]; % 75~100
rgb2=[1 0 0]; % 0~50
rgb3=[0 0 1]; % 50~75
rgb4=[1 1 0];

figure(3),
for i = 1 : Data_num
    hold on;
    if (practical_pulser_num50(i) == 0 && practical_pulser_num75(i) ==0 &&  practical_pulser_num100(i) ==0)
             plot(focus_hisx(i) , focus_hisy(i),'.' ,'Color' , [ rgb1(1)  rgb1(2)  rgb1(3)] ,'LineWidth',0.00005); %hold on;
    elseif (practical_pulser_num75(i) == 0 && practical_pulser_num100(i) ==0&& practical_pulser_num25(i) ==0)
             plot(focus_hisx(i) , focus_hisy(i),'.' ,'Color' , [ rgb2(1)  rgb2(2)  rgb2(3) ] ,'LineWidth',0.00005);% hold on;
    elseif (practical_pulser_num100(i) == 0 && practical_pulser_num25(i)==0&& practical_pulser_num50(i) ==0)
             plot(focus_hisx(i) , focus_hisy(i),'.' ,'Color' , [ rgb3(1)  rgb3(2)  rgb3(3) ] ,'LineWidth',0.00005); %hold on;
    elseif (practical_pulser_num25(i) == 0 && practical_pulser_num50(i)==0&& practical_pulser_num75(i) ==0)
             plot(focus_hisx(i) , focus_hisy(i),'.' ,'Color' , [ rgb4(1)  rgb4(2)  rgb4(3) ] ,'LineWidth',0.00005); %hold on;            
    end
    hold off;
end

for i = 1 : Data_num
    hold on;
    if (practical_pulser_num50(i) == 0 && practical_pulser_num75(i) ==0 &&  practical_pulser_num100(i) ==0)
             plot(-focus_hisx(i) , focus_hisy(i),'.' ,'Color' , [ rgb1(1)  rgb1(2)  rgb1(3)] ,'LineWidth',0.00005); %hold on;
    elseif (practical_pulser_num75(i) == 0 && practical_pulser_num100(i) ==0&& practical_pulser_num25(i) ==0)
             plot(-focus_hisx(i) , focus_hisy(i),'.' ,'Color' , [ rgb2(1)  rgb2(2)  rgb2(3) ] ,'LineWidth',0.00005);% hold on;
    elseif (practical_pulser_num100(i) == 0 && practical_pulser_num25(i)==0&& practical_pulser_num50(i) ==0)
             plot(-focus_hisx(i) , focus_hisy(i),'.' ,'Color' , [ rgb3(1)  rgb3(2)  rgb3(3) ] ,'LineWidth',0.00005); %hold on;
    elseif (practical_pulser_num25(i) == 0 && practical_pulser_num50(i)==0&& practical_pulser_num75(i) ==0)
             plot(-focus_hisx(i) , focus_hisy(i),'.' ,'Color' , [ rgb4(1)  rgb4(2)  rgb4(3) ] ,'LineWidth',0.00005); %hold on;
    end
    hold off            
end

for i = 1 : Data_num
    hold on;
    if (practical_pulser_num50(i) == 0 && practical_pulser_num75(i) ==0 &&  practical_pulser_num100(i) ==0)
             plot(-focus_hisx(i) , -focus_hisy(i),'.' ,'Color' , [ rgb1(1)  rgb1(2)  rgb1(3)] ,'LineWidth',0.00005); %hold on;
    elseif (practical_pulser_num75(i) == 0 && practical_pulser_num100(i) ==0&& practical_pulser_num25(i) ==0)
             plot(-focus_hisx(i) , -focus_hisy(i),'.' ,'Color' , [ rgb2(1)  rgb2(2)  rgb2(3) ] ,'LineWidth',0.00005);% hold on;
    elseif (practical_pulser_num100(i) == 0 && practical_pulser_num25(i)==0&& practical_pulser_num50(i) ==0)
             plot(-focus_hisx(i) , -focus_hisy(i),'.' ,'Color' , [ rgb3(1)  rgb3(2)  rgb3(3) ] ,'LineWidth',0.00005); %hold on;
    elseif (practical_pulser_num25(i) == 0 && practical_pulser_num50(i)==0&& practical_pulser_num75(i) ==0)
             plot(-focus_hisx(i) , -focus_hisy(i),'.' ,'Color' , [ rgb4(1)  rgb4(2)  rgb4(3) ] ,'LineWidth',0.00005); %hold on;
    end         
    hold off            
end

for i = 1 : Data_num
    hold on;
    if (practical_pulser_num50(i) == 0 && practical_pulser_num75(i) ==0 &&  practical_pulser_num100(i) ==0)
             plot(focus_hisx(i) , -focus_hisy(i),'.' ,'Color' , [ rgb1(1)  rgb1(2)  rgb1(3)] ,'LineWidth',0.00005); %hold on;
    elseif (practical_pulser_num75(i) == 0 && practical_pulser_num100(i) ==0&& practical_pulser_num25(i) ==0)
             plot(focus_hisx(i) , -focus_hisy(i),'.' ,'Color' , [ rgb2(1)  rgb2(2)  rgb2(3) ] ,'LineWidth',0.00005);% hold on;
    elseif (practical_pulser_num100(i) == 0 && practical_pulser_num25(i)==0&& practical_pulser_num50(i) ==0)
             plot(focus_hisx(i) , -focus_hisy(i),'.' ,'Color' , [ rgb3(1)  rgb3(2)  rgb3(3) ] ,'LineWidth',0.00005); %hold on;
    elseif (practical_pulser_num25(i) == 0 && practical_pulser_num50(i)==0&& practical_pulser_num75(i) ==0)
             plot(focus_hisx(i) , -focus_hisy(i),'.' ,'Color' , [ rgb4(1)  rgb4(2)  rgb4(3) ] ,'LineWidth',0.00005); %hold on;
    end         
    hold off            
end
axis square
axis image
%%
  
    title(' The number of pulser at each focus '); xlabel('lateral axis [mm]'); ylabel('elevational axis [mm]')
figure(4)  
for i = 1 : Data_num
    hold on;
             plot(focus_hisx(i) , focus_hisy(i), 'square' ,'Color' , [ practical_pulser_num01(i) practical_pulser_num01(i) practical_pulser_num01(i)] ,'LineWidth',4.3); %hold on;
    hold off           
end
for i = 1 : Data_num
    hold on;
             plot(focus_hisx(i) , -focus_hisy(i), 'square' ,'Color' , [ practical_pulser_num01(i) practical_pulser_num01(i) practical_pulser_num01(i)] ,'LineWidth',4.3); %hold on;
    hold off           
end
for i = 1 : Data_num
    hold on;
             plot(-focus_hisx(i) , focus_hisy(i), 'square' ,'Color' , [ practical_pulser_num01(i) practical_pulser_num01(i) practical_pulser_num01(i)] ,'LineWidth',4.3); %hold on;
    hold off           
end
for i = 1 : Data_num
    hold on;
             plot(-focus_hisx(i) , -focus_hisy(i),'square' ,'Color' , [ practical_pulser_num01(i) practical_pulser_num01(i) practical_pulser_num01(i)] ,'LineWidth',4.3); %hold on;
    hold off           
end    
    title(' The number of pulser at each focus '); xlabel('lateral axis [mm]'); ylabel('elevational axis [mm]')
xlim([-focus_hisx(xy_range*z_range)*1.00 focus_hisx(xy_range*z_range)*1.00] )
ylim([-focus_hisx(xy_range*z_range)*1.00 focus_hisx(xy_range*z_range)*1.00] )
axis square
colormap(gray); colorbar;  clim([0 64]);

%-------------------------------therapeutic quantization--------------------
p = 1;
delay_map = zeros(20+xy_range,20+z_range);
for k = 1:xy_range
    for i = 1:xy_range
    delay_map(41-k,10+i) = practical_pulser_num01(p);
    p = p + 1;
    end

end

for x= 1 : xy_range
    for y = 1:xy_range
        for xran = 1 : 10
            for yran = 1 : 10

                if delay_map(41-k,10+i) > delay_map(41-k-(5-xran),10+i-(5-yran))
                    delay_map(41-k,10+i) = delay_map(41-k-(5-xran),10+i-(5-yran));
                end
            end
        end
    end
end

delay_map30 = zeros(xy_range,z_range);
for x = 1:xy_range
    for y =1:xy_range
        delay_map30(x,y) = delay_map(10+x,10+y);
    end
end


%-------------------------------therapeutic quantization--------------------


for i = 1 : Data_num
    for k = 1 : 63
        delay_origin_diff(k ,i) = delay_origin(k,i) - delay_origin(k+1, i);
    end
end
%%% Pulser numbers



%%% MAX14808 clk frequency = 160MHz -> 6.25*1e-9
for i = 1:Data_num
    practical_num         = find ( delay_origin_diff(:,i) < 6.5*1e-9 );
    [practical_remove(1,i) , ~] = size(practical_num);
end


practical_remove01 = normalize(practical_remove,'range');
%%% MAX14808 clk frequency = 160MHz -> 6.25*1e-9


% figure(2) , 
% hold on

% for i = 1 : Data_num
% 
%     plot(focus_hisx(i) , focus_hisy(i),'.' ,'Color' , [ 0  practical_pulser_num01(i)  0 ] ...
%         ,'LineWidth',0.00005)
% 
% 
% end


%%% Numbers of Pulsers about each z-degree , graph
for i = 1:z_range

    Pulser_deg(i,:) = Pulser(xy_range*(i-1)+1 : xy_range*i);

end 
% 
% figure(2),
% for idx = 1: xy_range
% subplot(xy_range,z_range,idx)
% plot(1:xy_range , practical_pulser_num_cell(idx,:));
% ylim([min(Pulser) max(Pulser)])
% title([' # of Pulsers latitude =  ' num2str(idx), '[Deg]']);
% xlabel('longitude [Deg]')
% xticks([0:xy_range/2:xy_range]); xlim([0 xy_range]);
% ylabel('# of Pulsers')
% end
% %%% Numbers of Pulsers about each z-degree , graph







