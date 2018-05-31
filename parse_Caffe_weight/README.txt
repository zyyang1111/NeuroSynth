This file contains codes that parse the weight value output file from Caffe and save the weight information of each layer to CSV files sperately. 
1. The weight output is saved in a TXT file. You need to specifiy the path of this file. A sample output file (new_res.txt) is shown in this folder. The format is described as follows: 
	[#output feature maps, #input feature maps, dimesion of a filter]
        weight value in fixed point, weight value in integer format

2. If you would like to round round the weight value, please set "quantitize" (quantitize = 1). This will round the weight according to the property of value-driven optimized multipliers. In order to achieve this, you need to first generate a weight rouding file with the script in the "study_weights" folder. 