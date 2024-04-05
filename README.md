# AFM image simulation

This code is developed for simulating AFM image from atomic models. *Shall we upload the AFM data?*

Data: 03/2024.

This code is tested using Matlab 2023b and Python.  

For more information about the project, algorithms, and related publications please refer to the [Chen Group website](https://chenlab.matse.illinois.edu/). 

## Reference
If you find our codes are helpful to your publication, please cite:

Xiaoxu Li, Qing Guo, Yatong Zhao, Chang, Qian, Maxime Pouvreau, Trenton R Graham, Ping Chen1, Lili Liu, Chang Liu, Benjamin A Legg, Qian Chen, Aijun Miao, Zheming Wang, James J. De Yoreo, Carolyn I Pearce, Aurora E. Clark, Kevin M. Rosso, and Xin Zhang, "Dissolution of Mineral via Polynuclear Complexes" submitted to _Science_

## Getting started

### 1. Generate dataset
This part of code is developed to simulate AFM images of randomly generated atomic model of gibbsite (001) surface with defects and adatoms. The codes can be used upon downloading, data of 'Adatom libirary.mat' and 'gibbsite.vasp' are included.
- 1.1. Run 'AtomicSimulation_perfect.m' to obtain extended atomic model of gibbsite (001) surface as 'output/Perfect structure/Perfect structres_x.mat'. The code can generate a series of structures of different in plane rotation.
- 1.2. Run 'AtomicSimulation_defect.m' to obtain atomic models with defects and adatoms as 'output/Defected structure/Defected structures_x-y.mat'. For each perfect structure, the code can generate multiple random defect structures.
- 1.3. Run 'ForceCalculation_calcMatrix.m'. This step is generate matrix with all the surface information, and accelerate the next step of force calculation with randomized tip parameters. The code can go through all the defected structures in 'output/Defected structure/' and save output in 'output/Forces matrix/'.
- 1.4. Run 'DatasetGeneration.m' for actual dataset generation, or run 'DatasetGeneration_demo.m' to generate a few examples and learn about the workflow.

### 2. U-Net training and prediction
This part of code is developed to train a U-Net for predicting Al atom position of gibbsite (001) surface. The generated dataset from the last step ('1. Generate dataset/output/Dataset/image_stack.h5' can be directly used for the training, which is renamed as 'Training set+validation set.h5' as uploaded.
- 2.1 Install python environment a listed in 'tf gpu env.txt'.
	- conda create -n tf-gpu python=3.9.7
	- conda activate tf-gpu
	- conda install tensorflow-gpu=2.6.0
	- conda install keras
	- pip install matplotlib
- 2.2 Run 'Training U-Net.ipynb' if you want to train the U-Net model. The prediction of validation set will be generated from this file.
- 2.3 A pretrained U-Net model is also uploaded, by running 'Inspection U-Net.ipynb' to test the performance. You need some real AFM images with a input shape of (n,128,128) to run this.

### 3. Predict Al position
This part of code is developed to get Al atom position by upsampling. The prediction from U-Net as 'prediction_validation_20240312.h5' can be directly used.
- 3.1 Run 'Prediction2AlPosition_demo.m' to learn about the upsampling patterns, workflow, and the tracking error estimation.
- 3.2 To get all Al positions for real datasets, you need 'prediction_xxx.h5', and run 'Prediction2AlPosition_dataset.m'.




