# Integrated System for Seismic Response Monitoring and Risk Assessment of Buildings

## Project Overview
 
This project develops a comprehensive ****Integrated System for Seismic Response Monitoring and Risk Assessment of Buildings** (ISRAB)**. The system integrates real-time sensor data with deep learning models (Improved Deep Embedded Clustering, IDEC) to effectively predict changes in building states and assess risk. By incorporating **signal transmission and early warning**, **impact estimation** and **earthquake risk assessment**, the system provides full-cycle early warnings, evaluations, and decision support for earthquake disasters.

The project consists of three main modules:
1. **Signal Transmission and Early Warning**: Real-time sensor data acquisition, prediction, and alarm mechanisms.
2. **Impact Estimation**: Using building damage simulations, and multithreaded computations to evaluate potential building damage and fatalities.
3. **Earthquake Risk Assessment**: In conjunction with Impact Estimation, estimating potential earthquake event times based on fault data and the Gutenberg-Richter (GR) Law to provide foundational information for risk assessment. 


## Code Functionality Detailed Explanation

### 1. **signal_transmission _and_early_warning**

This folder handles real-time data transmission, storage, prediction, and alarm functionalities, mainly processing sensor data acquisition and damage prediction.

#### **`signal_decode/read_serial.py`**
- Responsible for reading sensor acceleration data from the serial port, converting it into a float array. The serial port address and other communication parameters must be configured according to the user's sensor setup and network configuration. These serial port settings can be customized in the **config** for easier modification. Example data is also provided in the MySQL database table for reference.
- Checks the data’s peak value and returns the corresponding alarm level through **`check_alarm_level()`**.

#### **`signal_decode/alarm.py`**
- Checks the alarm level based on the peak value of the acceleration data and provides four alarm levels: `Normal`, `Caution`, `Warning`.

#### **`upload/mysql_connector.py`**
- Inserts real-time sensor data into the **MySQL** database for storage, providing data support for subsequent analysis. Example data is provided in the MySQL database table for reference.


#### **`early_warning/prediction.py`**
- Loads deep learning models from the **model** folder and predicts building damage based on acceleration data. The **model** folder also contains the IDEC training process and network structure used for the prediction.
- Returns prediction results and the alarm level based on the prediction’s result.


#### **`config/config.py`**
- Stores configuration information, such as serial port settings, database configuration, and alarm thresholds.
- Users can modify this file to match different hardware or system requirements.

### 2. **impact_estimation**

This folder is responsible for impact estimation, simulating the damage to buildings during an event, and evaluating potential building damage, fatalities, and repair time.

#### **`Example_Estimation.py`**
- Provides an example for impact estimation, showing how to evaluate building damage and fatalities based on seismic intensity and damage models.

#### **`Gme.py`**
- Responsible for estimating seismic intensity and providing the necessary data for building damage assessment within the scope of risk analysis.
- Evaluates the seismic effects on buildings through seismic intensity models, contributing to the overall risk analysis process.

#### **`Lsc.py`**
- Uses **Gaussian Copula** and **Monte Carlo sampling** techniques to generate quantiles of correlated random variables related to building damage.
- The **Kriging interpolation** method is applied to improve computational efficiency.

#### **`Rscrv.py`**
- Uses **multithreading** to inverse-transform the sampled random variables and determine potential building loss and fatalities based on risk exposure data.

### 3. **earthquake_risk**

This folder focuses on earthquake risk assessment, particularly the prediction of potential earthquake times based on faults and the GR Law, providing foundational information for risk assessment.

#### **`gem_active_faults_harmonized.shp`**
- A **Shapefile** containing geographic data on potential faults. It provides essential information on seismic activity, including fault length, annual occurrence rate, and more. This data is global in scope, and users can select specific regions for targeted analysis, which is used in the risk assessment and earthquake event simulation.

#### **`potential_earthquake_event_stimulation.py`**
- Simulates potential earthquake event times based on the fault data and the **Gutenberg-Richter (GR) Law**, providing foundational data for the risk assessment process. Using earthquake event records from the USGS and the GR law, it determines the probability of different earthquake magnitudes, providing essential data for risk assessment.

## Frontend and API Integration

### 1. **databsase**

This part handles data storage and provides the necessary database examples for the ISRAB system.

- The folder **mysql_table** contains database examples relevant to ISRAB. These databases need to be customized by the user; however, for clarity and convenience in demonstrating the system's functionality, example databases have been uploaded for reference.

### 2. **api configuration**

This part provides the interfaces for interacting with external systems through API calls.

- The folder **api** contains a JAR file, which includes all the APIs that handle the data transfer process from MySQL to the frontend page. If users wish to maintain consistency between the data structure and the frontend, they can directly configure and migrate the setup.


### 3. **frontend configuration**

This part showcases the application of a 3D-printed farmhouse in the ISRAB system, demonstrating the complete process of structural monitoring, early warning, consequence assessment, and risk analysis.

- The frontend pages are written in JavaScript. Users can directly reference these pages for deployment or modify and optimize them according to their own needs.

Framework of Integrated System for Seismic Response Monitoring and Risk Assessment of Buildings

### 4 . **figure**

This part showcases the application of a 3D-printed farmhouse in the ISRAB system, demonstrating the complete process of structural monitoring, early warning, consequence assessment, and risk analysis.