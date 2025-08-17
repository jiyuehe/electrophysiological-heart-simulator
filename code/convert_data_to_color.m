function [hue,map_color] = convert_data_to_color(data_max,data_min,data_threshold,data,color_wheel_type)

debug_plot = 0; % show the color scale
if debug_plot == 1
    hue = (1:360) / 360;
    
    h_s_v = [];
    h_s_v(:,1) = hue;
    h_s_v(:,2) = 1;
    h_s_v(:,3) = 1;
    map_color = hsv2rgb(h_s_v);
    
    image = zeros(360,200,3);
    for i = 1:size(image,2)
        image(:,i,1) = map_color(:,1);
        image(:,i,2) = map_color(:,2);
        image(:,i,3) = map_color(:,3);
    end
    
    image(1:10,:,:) = 0; % make the first few lines black
    
    figure;
    imshow(image);
end

data(data<data_min) = data_min;
data(data>data_max) = data_max;

if strcmp(color_wheel_type,'red_megenta')
    hue = (data - data_min) / (data_max - data_min) * 300/360;
elseif strcmp(color_wheel_type,'red_red')
    hue = (data - data_min) / (data_max - data_min) * 360/360;
elseif strcmp(color_wheel_type,'red_blue')
    hue = (data - data_min) / (data_max - data_min) * 240/360;
end

h_s_v = [];
h_s_v(:,1) = hue;
h_s_v(:,2) = 1;
h_s_v(:,3) = 1;
map_color = hsv2rgb(h_s_v);

% assign non-active ones gray color
non_active_id = data <= data_threshold;
map_color(non_active_id,:) = 0.8;

map_color(isnan(map_color)) = 0.8;

end