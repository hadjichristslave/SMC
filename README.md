# SMC

Sequential monte carlo method for dependent Dirichlet process.
The sampler is based on the paper as found in http://jmlr.org/proceedings/papers/v33/neiswanger14.pdf.

The model is similar to that as the sampler is created to cluster point clouds rather than pixels.

#Random forest

The model output from the sampler is given as input to a random forest decision tree.
The outcome of this tree is what the algorithm perceives as what we are currently seing through our kinect sensors.

#Dependencies

PCL point cloud operations library

Eigen

OpenCV

Libconfig

Sqlite3(3 is important)

The project It is compiled using -std=c++11 flag

#How to use

Main.cpp expects input in a file of format:

*x,y,z,Kullback-leibler,EMD,Hellinger,Bin0-Bin26*

X,Y,Z represent the positional information

KL,EMD,Hellinger the distance distribution between a point and its K nearest neighbors.

The colour spectrum is discretized in bins and colour counts of the neighbor pixels are passed as input.

For more information on the theory behind the sampler check [this](http://people.ee.duke.edu/~lcarin/Yunchen8.15.2014.pdf presentation) and the related paper.

A random forest decision layer is added to classify if every cluster is an object the model has encountered before.

That way matching in the distribution space is done. The main motivation is to use this layer in the landmark detection phase of SLAM methods. It will greatly reduce the dimensionality of point clouds as it introduces an very extensive reprentation of a point cloud.

//Todo, dynamic number of colour bins
The file path is specified in the config file

Linker options are

-lconfig++

-lgsl

-lgslcblas

-lopencv_core

-lopencv_highgui

-lopencv_imgproc

-lsqlite3

#External packages used

**sqlite3pp** *boost free version (source found @ github: https://github.com/hadjichristslave/sqlite3pp)*

**random forest** *by bjorn andres( Source found @ github: https://github.com/bjoern-andres/random-forest)*


#Output

The project outputs in file specified in config file the clusters found on every sample of every pixel.

#Benchmark

The sampler is quite fast.
For 11K points in the cloud and a configuration of 20 samples and 3 particles it takes ~3s to output landmark ID's @ a Intel(R) Core(TM) i5-3210M CPU @ 2.50GHz
As is expected the weight calculations are quite expensive but they can be optimized even further.


#Attention

The files are still modified daily, so any weird behaviour might be due to development issues.

#Results

If the input is a point cloud the resulting clustering will look like this:
![alt tag](images/initialCluster.png)

![](https://github.com/github/hadjichristslave/SMC/images/initialCluster.png)

The cluster concentration on the part of the point cloud that represents the chair is noticeable.



