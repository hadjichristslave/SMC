# SMC

Sequential monte carlo method for dependent Dirichlet process.
The sampler is based on the paper as found in http://jmlr.org/proceedings/papers/v33/neiswanger14.pdf.

The model is similar to that as the sampler is created to cluster point clouds rather than pixels.

#Dependencies

Eigen
OpenCV
Libconfig

The project It is compiled using -std=c++11 flag

#How to use

Main.cpp expects input in a file of format:

x,y,z,Kullback-leibler,EMD,Hellinger,Bin0,Bin1,Bin2,Bin3,Bin4,Bin5,Bin6,Bin7,Bin8,Bin9,Bin10,Bin11,Bin12,Bin13,Bin14,Bin15,Bin16,Bin17,Bin18,Bin19,Bin20,Bin21,Bin22,Bin23,Bin24,Bin25,Bin26

X,Y,Z represent the positional information
KL,EMD,Hellinger the distance distribution between a point and its K nearest neighbors.

The colour spectrum is discretized in bins and colour counts of the neighbor pixels are passed as input.

//Todo, dynamic number of colour bins
The file path is specified in the config file

Linker options are

 --lconfig++
 --lgsl
 --lgslcblas
 --lopencv_core
 --lopencv_highgui
 --lopencv_imgproc



#Output

The project outputs in file specified in config file the clusters found on every sample of every pixel.

#Attention

The files are still modified daily, so any weird behaviour might be due to development issues.


