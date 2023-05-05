CXX = g++
MPICXX = mpicxx
CXXFLAGS = -O2 -std=c++11
LIBS = -larmadillo -lm
INCS = -I ./common \
-I ./methods \
-I ./methods/DIRK \
-I ./methods/DIRK/instances \
-I ./methods/MRI \
-I ./methods/MRI/instances \
-I ./methods/legacy \
-I ./problems \
-I ./controllers \
-I ./drivers/AdaptiveStep \
-I ./drivers/PaperRuns 

GenericAdaptiveDriver.exe : ./drivers/AdaptiveStep/GenericAdaptiveDriver.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(INCS)
	
GenericPaperRunsDriver.exe : ./drivers/PaperRuns/GenericPaperRunsDriver.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(INCS)
	
clean :
	\rm -rf *.txt *.exe *~ *.dSYM ./output/*/*.csv ./postprocessing/output/*/plots/*.png ./postprocessing/output/*/data/*.csv 
