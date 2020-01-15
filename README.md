# Data assimilation using the ensemble optimal interpolation 

This is the code used in the published paper 
[(link)](https://pdf.sciencedirectassets.com/280203/1-s2.0-S1877050917X00185/1-s2.0-S1877050917324018/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEOv%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLWVhc3QtMSJHMEUCIF3HyEj4mBd52IK3zURu9duVe8bGG307y0ilSWDBO2ddAiEA8YWMCaHkavyyyVaCcua6nmxsMZWq7JQsa%2Fw7AU2c6G8qtAMIcxACGgwwNTkwMDM1NDY4NjUiDDnfD9In0KNUouQGkiqRA36VqpAJZXvfNnC2ooX6B8PC%2FVcZhkS82EMgsfgRTYEGqrZGqMJ%2BUKA3%2FlkbeJsoixCE1cl%2FijXhEJxehplFadwbJylGZm5Nt9K%2BFxRUYT%2BCEjv5MYX8Wz0mPvB4WSFERzd4TnN0wsXztgdlOaMMeM8m2RhvwyJa70Fg5irDeClgONDuHV6bjOtSy0WIhp8bywcijpx%2B5mVOs81FIAggaDlGS1EG27dlzncYLED8TMQqCp6cdii%2FRe66e%2FevpKgutDkG3ebvk4AabCkDDRlBgwTJXvCeB%2FzJAY%2Bua8mXWKerqnlqCEGob%2FJChxHY9aC3ZvbW3nvA7MXc3rujoNmxdlkjjaylk88QzwRzBvrvSjaZTntCCEZDh9Of9%2Fpq7lr2cfCkqFTDGYmOZYMqfplfQLD8jrTzjuuugWrUJvy3ou01SOFJ8EmzDPRHEItmIXPyg9PPXheAhGTT7FRyw5DtrX01cnF3whlPsbH1kQhW1mKvxpfBAVnm9QvPH0FihpnkrbRgLHCPmXQQQCvjySsX1E%2FvMMTN%2B%2FAFOusBmmVohI%2FZF8R416EzBvOm6g0TeWqMgrDSIP%2F1iDG031qdPQ%2B58sEPO4ttDvQO7YXO8bLEQFPcxFvB2CFWcHKNYDLygugpeaNRQARoR2naqL61BNI8EPS7%2FM4XXPJ3pDV39xnVyUwILuOXv2cGPZRZj3nPDB8lKRet48FweeRROVHAPwT%2Bg7UNvYhsUtKzqWcFi2YITM1SpBbzDGphk%2BfiM3x1k4tPmUfQVKKdiCX%2BIjuQ9WEgVFWXfTnm5Q3UntVB1%2Br7uBeGP9%2FahDQx03M9UAxm9vOLTSHAhIpcF1tTjv%2BJx6LzINzcAxblCw%3D%3D&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20200115T112738Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTYSTT3PKHB%2F20200115%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=11104111acd2c5b5bbde5995b14ea96773c4a68f5d1420f133a2769de3d96376&hash=ef367e610d7f5354ba7e1b7e81a829c78afed3f7f5b528b52e75465ee62a5945&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=S1877050917324018&tid=spdf-5612e98e-b063-4c91-a35b-d67c416c63ce&sid=c6eb608775fd5340483b0a5-04e8b0c018afgxrqb&type=client)
and implemented in MATLAB.

## 1. Introduction
 **Ensemble optimal interpolation (EnOI)** is an efficient, having a relatively low computational cost yet powerful
 method for correcting numerical simulation models in according to the field measurements. Basic EnOI is used in different areas 
 of geoscience including meteorology and oceanography. 
 
 The goal was to apply the method of EnOI to post-process the numerically simulated wind field in the way that it would maximally 
 fit the actual measurements at the weather stations. The basic principle behind this technique is to interpolate discrepancies 
 between measurements and simulation following the spatiotemporal consistency 
 of the modeled wind field. However, its statistical properties are in general non-stationary, spatially inhomogeneous 
 as well as the quality of measurements vary from station to station. 
 Thus, such a correction can introduce undesirable artifacts (outliers) into the output. We adressed this problem by 
 adding the ridge regression into the basic EnOI scheme.
 
 The algorithm was applied to assimilate observations to the wind field simulated for the Polar Arctic region.
 

 
 ## 2. Mathematical model
 
### Ensemble optimal interpolation

Resulting correction ![equation](https://latex.codecogs.com/gif.latex?x_i%5E%7Ba%7D) can be expressed as a weighted linear combination 
of the differences between observed and modeled time series 
![equation](https://latex.codecogs.com/png.latex?%5CDelta_n%28%5Ctau%29%20%3D%20y_n%5E%7Bo%7D%28%5Ctau%29-Hx_n%5E%7Bb%7D%28%5Ctau%29) 
added then to each grid point of the modeled (background) field ![equation](https://latex.codecogs.com/png.latex?x_i%5E%7Bb%7D):

![equation](https://latex.codecogs.com/png.latex?x_i%5E%7Ba%7D%28t%29%20%3D%20x_i%5E%7Bb%7D%28t%29%20&plus;%20%5Csum_%7Bn%3D1%7D%5E%7BN%7D%20%5Csum_%7B%5Ctau%3Dt-T%7D%5E%7Bt&plus;T%7D%20w_n%5E%7B%28i%2Ct%29%7D%28%5Ctau%29%5CDelta_n%28%5Ctau%29)

where:\
![equation](https://latex.codecogs.com/png.latex?N) - number of observation points,\
![equation](https://latex.codecogs.com/png.latex?t) - time instance of the background model,\
![equation](https://latex.codecogs.com/png.latex?T) - maximum time lag that depends on the variability of the studied process.

The conventional OI results in the Best Linear Unbiased Estimate (BLUE) ![equation](https://latex.codecogs.com/png.latex?%5Ctextbf%7Bx%7D%5E%7Ba%7D) of the true system state column vector 
![equation](https://latex.codecogs.com/png.latex?%5Ctextbf%7Bx%7D%5E%7Bt%7D) having available vector of observations ![equation](https://latex.codecogs.com/png.latex?%5Ctextbf%7By%7D%5E%7Bo%7D), 
initial background field state vector ![equation](https://latex.codecogs.com/png.latex?%5Ctextbf%7Bx%7D%5E%7Bb%7D), and estimated error covariances:

![equation](https://latex.codecogs.com/png.latex?%5Ctextbf%7Bx%7D%5E%7Ba%7D%3D%5Ctextbf%7Bx%7D%5E%7Bb%7D&plus;%5Ctextbf%7BK%7D%28%5Ctextbf%7By%7D%5E%7Bo%7D-%5Ctextbf%7BH%7D%5Ctextbf%7Bx%7D%5E%7Bb%7D%29)

![equation](https://latex.codecogs.com/png.latex?%5Ctextbf%7BK%7D%3D%5Ctextbf%7BP%7D%5E%7Bb%7D%5Ctextbf%7BH%7D%5E%7BT%7D%28%5Ctextbf%7BH%7D%5Ctextbf%7BP%7D%5Eb%5Ctextbf%7BH%7D%5E%7BT%7D&plus;%5Ctextbf%7BR%7D%29%5E%7B-1%7D)

where:\
![equation](https://latex.codecogs.com/png.latex?%5Ctextbf%7BP%7D%5Eb%20%3D%20%5Clangle%20%28%5Ctextbf%7Bx%7D%5Eb%20-%20%5Ctextbf%7Bx%7D%5Et%29%7B%28%5Ctextbf%7Bx%7D%5Eb%20-%20%5Ctextbf%7Bx%7D%5Et%29%7D%5ET%20%5Crangle%24%20and%20%24%5Ctextbf%7BR%7D%20%3D%20%5Clangle%20%28%5Ctextbf%7By%7D%5Eo%20-%20%5Ctextbf%7Bx%7D%5Et%29%7B%28%5Ctextbf%7By%7D%5Eo%20-%20%5Ctextbf%7Bx%7D%5Et%29%7D%5ET%20%5Crangle) are the matrices of estimated covariances of the background field and observation errors, respectively;\
![equation](https://latex.codecogs.com/png.latex?%5Ctextbf%7BH%7D) - the representation matrix of mapping observations onto the corresponding model grid;\
![equation](https://latex.codecogs.com/png.latex?%5Cmathbf%7BK%7D) - the weight matrix referred to as the Kalman gains, since the OI corresponds to the analysis stage of the Kalman filter. 
