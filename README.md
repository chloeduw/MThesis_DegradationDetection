# MThesis_DegradationDetection
This repository contains the files used in the scope of the following Master's Thesis (Bruface): Degradation Detection and Localisation in Battery Packs
At the beginning of this work, the **EHM simulator files** were provided. The rest of the files were created in the scope of this thesis.  

---

## 📂 EHM Simulator

This folder contains the following files (as given at the beginning of this work):

- `script_simu_ehm.m`: Script applying the EHM simulator  
- `EHM_fcn_model.m`: Function applying the EHM simulator  
- `output_EHM.m`: Function applying Eq. (EHM_y)  
- `EHMlinapprox.m`: Function applying Eq. (EHMlin)  
- `init_cs.m`: Function computing the initial solid concentrations from the initial voltage and the parameters of the cell  
- `params_LCO.m`: Structure containing the parameters of the cell  
- `refPotentialCathode.m`: Function giving the reference potential of the cathode and its derivative in function of the CSC  
- `refPotentialAnode.m`: Function giving the reference potential of the anode and its derivative in function of the CSC  

---

## 📂 Degradation Detection Files

### 🔹 EKF  
- `ekf.m`: Function applying the EKF algorithm  
- `deta.m`: Function applying the derivative of the electrode overpotentials  
- `JacobianMeasMat.m`: Function calculating the Jacobian of Eq. (EHM_y)  
- Functions from the EHM simulator (some may have been modified)  

### 🔹 MMAE  
- `conditionalDensityFct.m`: Function applying Eq. (CondDstyFct)  

### 🔹 Simulink Block  
- `+EHMCell/`: Repository containing the `.ssc` file of the Simulink block and its icon  
- `EHMCell_lib.slx`: Library containing the Simulink block  
- `TestEHMCell.slx`: Simulink model (fresh 1-cell configuration)  
- `TestEHMCellDegraded.slx`: Simulink model (degraded 1-cell configuration)  
- `TestAgeingImpact_Intermediate_CC_Cycle.slx`: Simulink model (2 parallel cells: 1 fresh, 1 degraded)  
- `simu_effet_ageing.m`: Script used in data generation results  
- `EKFvalidation.m`: Script used in EKF validation results  
- `simu_ageing_intermediaire.m`: Script used in ageing detection results (two EKFs case)  
- `simu_ageing_intermediaire_3HYP.m`: Script used in ageing detection results (three EKFs case)  

---
