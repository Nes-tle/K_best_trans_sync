This code repository implements the following ICCV 2019 paper:
K-Best Transformation Synchronization. Yifan Sun, Jiacheng Zhuo, Arnav Mohan, Qixing Huang; The IEEE International Conference on Computer Vision (ICCV), 2019, pp. 10252-10261 
% Sample code:

>>Results = k_best_sync(pairmatches, Para);
>>save_pcs(pc, Results, 'scans/');

The code requires a minimum number of parameters. The code works for both dense graphs and sparse graphs. 

If you use the code for comparison, please cite the paper accordingly. 
@InProceedings{Sun_2019_ICCV,
author = {Sun, Yifan and Zhuo, Jiacheng and Mohan, Arnav and Huang, Qixing},
title = {K-Best Transformation Synchronization},
booktitle = {The IEEE International Conference on Computer Vision (ICCV)},
month = {October},
year = {2019}
} 
