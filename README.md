# GOI
Retrieving the Gradients of O-information (GOI) from data.
Current version is only for a zero-lag analysis, which works both for time series and for data collected across multiple subjects/trials; dynamical analysis will be implemented in the future (feel free to contribute).

It uses the gaussian normalization presented in Robin Ince's repository for gaussian copula mutual information (https://github.com/robince/gcmi)

The main files are 
- goi_gradients.m   computes exhaustively the gradients until a max order
- goi_oinfo.m      computes the O-information of the input 


### References
1. Scagliarini, Tomas, et al. "Gradients of O-information: Low-order descriptors of high-order dependencies." Physical Review Research 5.1 (2023): 013025.
2. Ince, Robin AA, et al. "A statistical framework for neuroimaging data analysis based on mutual information estimated via a gaussian copula." Human brain mapping 38.3 (2017): 1541-1573.
