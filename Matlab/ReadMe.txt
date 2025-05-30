==============================
README - MPEE_V1 (May 2025)
==============================

Author: Siyu Zhu
Institution: University of Oklahoma

------------------------------
Function: MPEE_V1.m
------------------------------
This function extracts 3D precipitation events from a gridded data cube using a morphological, seed-based method.

Input:
- PMinput: a 3D matrix (lat x lon x time) of hourly precipitation (in mm/h)

Optional parameters (with default values):
- 'SeedThreshold'       = 1.5    % High threshold for event core
- 'PeripheralThreshold' = 0.1    % Low threshold for peripheral growth
- 'FilterKernel'        = [5 5 5]% 3D smoothing filter size
- 'DilationIterations'  = 50     % Number of 3D dilation steps
- 'DilationRadius'      = 1      % Radius of spherical dilation
- 'MinEventSize'        = 200    % Minimum number of voxels to keep

Output:
- Eventlist: a struct array with fields:
    * id  - event ID
    * MS  - row indices (latitude)
    * NS  - column indices (longitude)
    * TS  - time indices (hours)
    * PS  - precipitation values (mm/h)

------------------------------
Demo usage:
------------------------------
1. Open `demo_main.m` and run it.
2. It loads test data `PMall_Vtest.mat` from ./TestData.
3. It calls MPEE_V1 and prints summary info for the first 5 events.

------------------------------
Notes:
------------------------------
- Make sure the folder `.\MPEE_V1` is added to the MATLAB path.
- This code is designed for hourly precipitation inputs.

For any questions, please contact: zhusiyu2023@gmail.com

