function [fig_size, file_array] = determine_NetInfo(str)

% Zhiyuan Yang   May 28 2017
% determine the figure size of each layer and the name of each layer 

if strcmp(str, 'AlexNet_caffe')
    fig_size = [55;
                27;
                13;
                13;
                13;
                1;
                1;
                1];
     file_array = {'conv1';
                   'conv2';
                   'conv3';
                   'conv4';
                   'conv5';
                   'fc6';
                   'fc7';
                   'fc8'};
elseif strcmp(str, 'GoogleNet_caffe')
    fig_size = [32;
                15*ones(2,1);
                7*ones(12,1);
                3*ones(30,1);
                1*ones(12,1);
                1];
    file_array = {'conv1_7x7_s2';
                  'conv2_3x3_reduce';
                  'conv2_3x3';
                  'inception_3a_1x1';
                  'inception_3a_3x3_reduce';
                  'inception_3a_3x3';
                  'inception_3a_5x5_reduce';
                  'inception_3a_5x5';
                  'inception_3a_pool_proj';
                  'inception_3b_1x1';
                  'inception_3b_3x3_reduce';
                  'inception_3b_3x3';
                  'inception_3b_5x5_reduce';
                  'inception_3b_5x5';
                  'inception_3b_pool_proj';
                  'inception_4a_1x1';
                  'inception_4a_3x3_reduce';
                  'inception_4a_3x3';
                  'inception_4a_5x5_reduce';
                  'inception_4a_5x5';
                  'inception_4a_pool_proj';
                  'inception_4b_1x1';
                  'inception_4b_3x3_reduce';
                  'inception_4b_3x3';
                  'inception_4b_5x5_reduce';
                  'inception_4b_5x5';
                  'inception_4b_pool_proj';
                  'inception_4c_1x1';
                  'inception_4c_3x3_reduce';
                  'inception_4c_3x3';
                  'inception_4c_5x5_reduce';
                  'inception_4c_5x5';
                  'inception_4c_pool_proj';
                  'inception_4d_1x1';
                  'inception_4d_3x3_reduce';
                  'inception_4d_3x3';
                  'inception_4d_5x5_reduce';
                  'inception_4d_5x5';
                  'inception_4d_pool_proj';
                  'inception_4e_1x1';
                  'inception_4e_3x3_reduce';
                  'inception_4e_3x3';
                  'inception_4e_5x5_reduce';
                  'inception_4e_5x5';
                  'inception_4e_pool_proj';
                  'inception_5a_1x1';
                  'inception_5a_3x3_reduce';
                  'inception_5a_3x3';
                  'inception_5a_5x5_reduce';
                  'inception_5a_5x5';
                  'inception_5a_pool_proj';
                  'inception_5b_1x1';
                  'inception_5b_3x3_reduce';
                  'inception_5b_3x3';
                  'inception_5b_5x5_reduce';
                  'inception_5b_5x5';
                  'inception_5b_pool_proj';
                  'loss3_classifier'};
else
    fprintf('Unknown Network!\n');
    fig_size = []; 
    file_array = [];
end