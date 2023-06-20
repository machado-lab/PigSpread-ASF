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

```
@article{sykes2023estimating,
  title={Estimating the effectiveness of control actions on African swine fever transmission in commercial swine populations in the United States},
  author={Sykes, Abagael L and Galvis, Jason A and O’Hara, Kathleen C and Corzo, Cesar and Machado, Gustavo},
  journal={Preventive Veterinary Medicine},
  pages={105962},
  year={2023},
}
```
**Abstract**
Given the proximity of African swine fever (ASF) to the U.S., there is an urgent need to better understand the possible dissemination pathways of the virus within the U.S. swine industry and to evaluate mitigation strategies. Here, we extended PigSpread, a farm-level spatially-explicit stochastic compartmental transmission model incorporating six transmission routes including between-farm swine movements, vehicle movements, and local spread, to model the dissemination of ASF. We then examined the effectiveness of control actions similar to the ASF national response plan. The average number of secondary infections during the first 60 days of the outbreak was 49 finisher farms, 17 nursery farms, 5 sow farms, and less than one farm in other production types. The between-farm movements of swine were the predominant route of ASF transmission with an average contribution of 71.1%, while local spread and movement of vehicles were less critical with average contributions of 14.6% and 14.4%. We demonstrated that the combination of quarantine, depopulation, movement restrictions, contact tracing, and enhanced surveillance, was the most effective mitigation strategy, resulting in an average reduction of 79.0% of secondary cases by day 140 of the outbreak. Implementing these control actions led to a median of 495,619 depopulated animals, 357,789 diagnostic tests, and 54,522 movement permits. Our results suggest that the successful elimination of an ASF outbreak is likely to require the deployment of all control actions listed in the ASF national response plan for more than 140 days, as well as estimating the resources needed for depopulation, testing, and movement permits under these controls.
