CXX = g++
CXXFLAGS = -O3 -std=c++17 -Wall
LDFLAGS = -lgmp

EIGEN_PATH = /usr/include/eigen3
CXXFLAGS += -I$(EIGEN_PATH)

# Targets
all: collatz ulam-gap correlation fractional-position

# Main Simulation (Paper 1 & 2)
collatz: collatz.cpp
	$(CXX) $(CXXFLAGS) -o collatz collatz.cpp $(LDFLAGS)

# Fractional position data (Paper 2)
fractional-position: fractional-position.cpp
	$(CXX) $(CXXFLAGS) -o fractional-position fractional-position.cpp $(LDFLAGS)

# Autocorrelation data (Paper 2)	
correlation: correlation.cpp
	$(CXX) $(CXXFLAGS) -o correlation correlation.cpp $(LDFLAGS)

# Spectral Gap / Ulam Method (Supplement)
ulam-gap: ulam-gap.cpp
	$(CXX) $(CXXFLAGS) -o ulam-gap ulam-gap.cpp $(LDFLAGS)

clean:
	rm -f collatz ulam-gap correlation fractional-position

ulam-gap-out.csv: ulam-gap
	./ulam-gap > ulam-gap-out.csv

autocorrelation.csv: correlation
	./correlation > autocorrelation.csv

fractional-position.csv: fractional-position
	./fractional-position > fractional-position.csv

# This will take while
generate-data: fractional-position.csv autocorrelation.csv ulam-gap.csv
	echo "Data: Done."

autocorrelation.eps: autocorrelation.csv
	Rscript plot-autocorrelation.R

benford-dist.eps: fractional-position.csv
	Rscript plot-fraction-position.R

spectral-gap.eps: ulam-gap-out.csv
	Rscript plot-spectral-gap.R

generate-plots: autocorrelation.eps benford-dist.eps spectral-gap.eps
	echo "Plots: Done."

clean-data:
	rm -f *.csv

clean-plots:
	rm -f *.eps
