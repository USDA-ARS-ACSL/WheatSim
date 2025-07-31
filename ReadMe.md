How to compile the "WheatSim" crop model:
1. Using the graphics card processor is optional (following CUDA steps below). This will speed up the model somewhat and is proof of concept.
2. Download "CUDA (Runtime) Toolkit" here "https://developer.nvidia.com/cuda-downloads". The CUDA version used for WheatSim development is 12.5. "CUDA (Runtime) Toolkit" is not "CuDNN", and "CuDNN" is not needed for WheatSim.
3. Install "CUDA (Runtime) Toolkit" and that will automatically install Nsight compiler and an API for your Visual Studio.
4. After you installed "CUDA (Runtime) Toolkit", open Visual Studio->new project, you will see that "cuda runtime" option. That means you are ready to go.
5. Download WheatSim crop model and compile.
6. There is an example using data from APSIM you can run from visual studio by modifying the debug options or run from the command line as: wheatsim .\ApsimExample\APS02\runAPS02TOS1HighN.dat
   
