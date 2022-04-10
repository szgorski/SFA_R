# Original description of the algorithm by Patrick Schäfer, Ulf Leser (2017) 
# in "Fast and Accurate Time Series Classification with WEASEL" 
# (Word ExtrAction for time SEries cLassification), 
# DOI: 10.1145/3132847.3132980, arXiv:1701.07681

# This implementation by Szymon Górski (2021 - 2022)

# Original implementation of WEASEL algorithm by Patrick Schäfer (2017 - 2021) 
# under GNU General Public License v3.0,
# https://github.com/patrickzib/SFA

# Based on Python implementation by Samuel Harford (2017 - 2018) 
# under GNU General Public License v3.0,
# https://github.com/patrickzib/SFA_Python

# Testing datasets from 2015 UCR Time Series Classification Archive
# available at: https://www.cs.ucr.edu/~eamonn/time_series_data/

# Testing datasets from 2018 UCR Time Series Classification Archive
# available at: https://www.cs.ucr.edu/~eamonn/time_series_data_2018/

source("WEASELClassifier.R")
options(max.print = 100000)

dataset2015Names <- c("50words", "Adiac", "ArrowHead", "Beef", "BeetleFly", 
                      "BirdChicken", "Car", "CBF", "ChlorineConcentration", 
                      "CinC_ECG_torso", "Coffee", "Computers", "Cricket_X", 
                      "Cricket_Y", "Cricket_Z", "DiatomSizeReduction", 
                      "DistalPhalanxOutlineAgeGroup", 
                      "DistalPhalanxOutlineCorrect", "DistalPhalanxTW", 
                      "Earthquakes", "ECG200", "ECG5000", "ECGFiveDays", 
                      "ElectricDevices", "FaceAll", "FaceFour", "FacesUCR", 
                      "FISH", "FordA", "FordB", "Gun_Point", "Ham", 
                      "HandOutlines", "Haptics", "Herring", "InlineSkate", 
                      "InsectWingbeatSound", "ItalyPowerDemand", 
                      "LargeKitchenAppliances", "Lighting2", "Lighting7", 
                      "MALLAT", "Meat", "MedicalImages", 
                      "MiddlePhalanxOutlineAgeGroup", 
                      "MiddlePhalanxOutlineCorrect", "MiddlePhalanxTW", 
                      "MoteStrain", "NonInvasiveFatalECG_Thorax1", 
                      "NonInvasiveFatalECG_Thorax2", "OliveOil", "OSULeaf", 
                      "PhalangesOutlinesCorrect", "Phoneme", "Plane", 
                      "ProximalPhalanxOutlineAgeGroup", 
                      "ProximalPhalanxOutlineCorrect", "ProximalPhalanxTW", 
                      "RefrigerationDevices", "ScreenType", "ShapeletSim", 
                      "ShapesAll", "SmallKitchenAppliances", 
                      "SonyAIBORobotSurface", "SonyAIBORobotSurfaceII", 
                      "StarLightCurves", "Strawberry", "SwedishLeaf", "Symbols",
                      "synthetic_control", "ToeSegmentation1", 
                      "ToeSegmentation2", "Trace", "TwoLeadECG", "Two_Patterns",
                      "UWaveGestureLibraryAll", "uWaveGestureLibrary_X", 
                      "uWaveGestureLibrary_Y", "uWaveGestureLibrary_Z", "wafer",
                      "Wine", "WordsSynonyms", "Worms", "WormsTwoClass", "yoga")
dataset2018Names <- c("ACSF1", "Adiac", "AllGestureWiimoteX", 
                      "AllGestureWiimoteY", "AllGestureWiimoteZ", "ArrowHead", 
                      "Beef", "BeetleFly", "BirdChicken", "BME", "Car", "CBF", 
                      "Chinatown", "ChlorineConcentration", "CinCECGTorso", 
                      "Coffee", "Computers", "CricketX", "CricketY", "CricketZ",
                      "Crop", "DiatomSizeReduction", 
                      "DistalPhalanxOutlineAgeGroup",
                      "DistalPhalanxOutlineCorrect", "DistalPhalanxTW", 
                      "DodgerLoopDay", "DodgerLoopGame", "DodgerLoopWeekend", 
                      "Earthquakes", "ECG200", "ECG5000", "ECGFiveDays", 
                      "ElectricDevices", "EOGHorizontalSignal", 
                      "EOGVerticalSignal", "EthanolLevel", "FaceAll", 
                      "FaceFour", "FacesUCR", "FiftyWords", "Fish", "FordA", 
                      "FordB", "FreezerRegularTrain", "FreezerSmallTrain", 
                      "Fungi", "GestureMidAirD1", "GestureMidAirD2", 
                      "GestureMidAirD3", "GesturePebbleZ1", "GesturePebbleZ2", 
                      "GunPoint", "GunPointAgeSpan", "GunPointMaleVersusFemale",
                      "GunPointOldVersusYoung", "Ham", "HandOutlines", 
                      "Haptics", "Herring", "HouseTwenty", "InlineSkate", 
                      "InsectEPGRegularTrain", "InsectEPGSmallTrain", 
                      "InsectWingbeatSound", "ItalyPowerDemand", 
                      "LargeKitchenAppliances", "Lightning2", "Lightning7", 
                      "Mallat", "Meat", "MedicalImages", "MelbournePedestrian", 
                      "MiddlePhalanxOutlineAgeGroup", 
                      "MiddlePhalanxOutlineCorrect", "MiddlePhalanxTW", 
                      "MixedShapesRegularTrain", "MixedShapesSmallTrain", 
                      "MoteStrain", "NonInvasiveFetalECGThorax1", 
                      "NonInvasiveFetalECGThorax2", "OliveOil", "OSULeaf", 
                      "PhalangesOutlinesCorrect", "Phoneme", 
                      "PickupGestureWiimoteZ", "PigAirwayPressure", 
                      "PigArtPressure", "PigCVP", "PLAID", "Plane", 
                      "PowerCons", "ProximalPhalanxOutlineAgeGroup", 
                      "ProximalPhalanxOutlineCorrect", "ProximalPhalanxTW", 
                      "RefrigerationDevices", "Rock", "ScreenType", 
                      "SemgHandGenderCh2", "SemgHandMovementCh2", 
                      "SemgHandSubjectCh2", "ShakeGestureWiimoteZ", 
                      "ShapeletSim", "ShapesAll", "SmallKitchenAppliances", 
                      "SmoothSubspace", "SonyAIBORobotSurface1", 
                      "SonyAIBORobotSurface2", "StarLightCurves", "Strawberry", 
                      "SwedishLeaf", "Symbols", "SyntheticControl", 
                      "ToeSegmentation1", "ToeSegmentation2", "Trace", 
                      "TwoLeadECG", "TwoPatterns", "UMD", 
                      "UWaveGestureLibraryAll", "UWaveGestureLibraryX", 
                      "UWaveGestureLibraryY", "UWaveGestureLibraryZ", 
                      "Wafer", "Wine", "WordSynonyms", "Worms", 
                      "WormsTwoClass", "Yoga")

catch <- list()
catch <- uvLoad(dataset2015Names[[4]])
train <- catch[[1]]
test <- catch[[2]]
classifier <- createClassifier()
classifier <- evaluate(classifier, train, test)
