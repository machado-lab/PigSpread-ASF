# PigSpread-ASF

Here we have developed a stochastic farm-level compartmental transmission model to simulate the spread of African swine fever (ASF) within the United States (U.S.). The model takes into account real population, animal movement and vehicle movement data to allow us to incorporate 6 transmission routes:
* Between-farm swine movements
* Vehicles moving swine from farm to farm
* Vehicles moving swine from farms to slaughterhouses
* Vehicles moving crew from farm to farm
* Vehicles delivering feed to farms
* Local spread (transmission via geogrpahic vicinity) 

As part of the model we assessed the effectiveness of 5 control actions for ASF:
1. Depopulation of detected farms 
2. 72 hour movement standstill
3. Contact tracing of direct (swine movements) and indirect (vehicle movements) contacts
4. Depopulation of direct contacts
5. Implementation of control areas and surveillance zones, including diagnostic testing and movement permits. 

In this repository, you will find the generic code to run the model, including calibration of parameters and implementation of the controls. 
