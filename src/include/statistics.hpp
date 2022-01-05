// this is a header for aggregate statistics,
// for mean, std, correlations
// that do not store all the data points
// useful if we have too much data and little memory


#ifndef STATS_HPP
#define STATS_HPP

#include <vector>
#include <string>
#include <cmath>


// most algorithms are from
// 
// Optimization of Pearson correlation coefficient
// calculation for DPA and comparison of different
// approaches
// 
// Petr Socha, Vojtěch Miškovský, Hana Kubátová, Martin Novotný


// mean

template <typename T>
class StatAggregatorMean {

public:
	StatAggregatorMean()
		: Mean(0), nSamples(0) {}

	virtual ~StatAggregatorMean() {}


	void ConsumeDataPoint(const T &x) {
		++nSamples;
		Mean = Mean + (x - Mean) / nSamples;
	}

	void Reset() {
		Mean = 0;
		nSamples = 0;
	}

	T CurrentValue() const { return Mean; }
	unsigned CurrentNSamples const { return nSamples; }

private:
	T Mean;
	unsigned nSamples;
};

// standard deviation
// using Welford's one pass algorithm
// https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Welford's_online_algorithm

template <typename T>
class StatAggregatorWelfordStD {

public:
	StatAggregatorWelfordStD()
		: S(0), Mean(), nSamples(0) {}

	virtual ~StatAggregatorWelfordStD() {}


	void ConsumeDataPoint(const T &x) {
		T prev_mean = Mean.CurrentValue();
		Mean.ConsumeDataPoint(x);
		T cur_mean = Mean.CurrentValue();
		S = S + (x - prev_mean) * (x - cur_mean);
	}

	void Reset() {
		S = 0;
		nSamples = 0;
		Mean.Reset();
	}

	T CurrentValue() const { return sqrt(S / Mean.CurrentNSamples()) ; }
	T CurrentMean() const { return Mean.CurrentValue() ; }

private:
	T S;
	StatAggregatorMean<T> Mean;
};

// correlation
// computed from covariance
// https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Online


#endif	//
