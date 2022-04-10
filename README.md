# Time Series Classification in R

This repository contains an R (4.0.1) implementation of symbolic time series representation (**SFA** - **Symbolic Fourier Approximation**) and a univariate (**WEASEL** - **Word ExtrAction for time SEries cLassification**) time series model for series data analytics, originally designed and implemented by Patrick Schäfer. <ins>Detailed description of algorithms is provided in his publication [1] and GitHub repository [2].</ins> The solution was also based on Samuel Harford’s Python implementation [3].

### Datasets

The **main.R** and **timeSeries.R** files are prepared for loading training and testing datasets in a form of 2015 and 2018 UCR Time Series Classification Archive. For testing the algorithm, the "Beef" dataset was provided in **datasets** folder and is loaded by default. The output of both this and Samuel Harford’s Python implementation were compared in the **Test.pdf** file.

All data must be provided in a form:

*label<sub>1</sub> value<sub>1,1</sub> value<sub>1,2</sub> … value<sub>1,n</sub>*<br />*label<sub>2</sub> value<sub>2,1</sub> value<sub>2,2</sub> … value<sub>2,n</sub>*<br />…<br />*label<sub>m</sub> value<sub>m,1</sub> value<sub>m,2</sub> … value<sub>m,n</sub>*

### TODO

The solution will be improved in the future to decrease the amount of data stored and possibly the computation time by redesigning the sequence of operations.

### Links

 1. Patrick Schäfer, Ulf Leser (2017), "Fast and Accurate Time Series Classification with WEASEL", <br />[DOI: 10.1145/3132847.3132980](https://doi.org/10.1145/3132847.3132980), [arXiv:1701.07681](https://arxiv.org/abs/1701.07681)

 2. Original implementation of the WEASEL algorithm by Patrick Schäfer (2017 - 2021) under GNU General Public License v3.0, <br />[https://github.com/patrickzib/SFA](https://github.com/patrickzib/SFA)

 3. Based on Python implementation by Samuel Harford (2017 - 2018) under GNU General Public License v3.0, <br />[https://github.com/patrickzib/SFA_Python](https://github.com/patrickzib/SFA_Python)

 4. Testing datasets from 2015 "The UCR Time Series Classification Archive": <br />Yanping Chen, Eamonn Keogh, Bing Hu, Nurjahan Begum, Anthony Bagnall, Abdullah Mueen and Gustavo Batista (2015), <br />[https://www.cs.ucr.edu/~eamonn/time_series_data/](https://www.cs.ucr.edu/~eamonn/time_series_data/)

 5. Testing datasets from 2018 "The UCR Time Series Classification Archive": <br />Hoang Anh Dau, Eamonn Keogh, Kaveh Kamgar, Chin-Chia Michael Yeh, Yan Zhu, Shaghayegh Gharghabi , Chotirat Ann Ratanamahatana, Yanping Chen, Bing Hu, Nurjahan Begum, Anthony Bagnall , Abdullah Mueen, Gustavo Batista and Hexagon-ML (2019), <br />[https://www.cs.ucr.edu/~eamonn/time_series_data_2018/](https://www.cs.ucr.edu/~eamonn/time_series_data_2018/)