
#include <string>
#include <algorithm>
#ifdef _WIN32
#include <sys/timeb.h>
#include <sys/types.h>
#else //!_WIN32
#include <limits.h>
#include <sys/time.h>
#include <unistd.h>
#endif


namespace stats {
// Checks that the sequence is sorted, returns true if it's sorted, false
// otherwise
template<typename ForwardIt>
bool is_sorted(ForwardIt first, ForwardIt last){
  ForwardIt next = first;
  ++next;
  for (; next != last; ++next, ++first){
    if (*first > *next){
      return false;
    }
  }
  return true;
}

// Get the value representing the `percent` percentile of the
// sorted samples using linear interpolation
double PercentileValue(const std::vector<double> &samples, const double percent){
  assert(!samples.empty());
  if (samples.size() == 1){
    return samples.front();
  }
  assert(percent >= 0.0);
  assert(percent <= 100.0);
  assert(stats::is_sorted(samples.begin(), samples.end()));
  if (percent == 100.0){
    return samples.back();
  }
  // Find the two nearest percentile values and linearly
  // interpolate between them
  const double rank = percent / 100.0 * (static_cast<double>(samples.size()) - 1.0);
  const double low_rank = std::floor(rank);
  const double dist = rank - low_rank;
  const size_t k = static_cast<size_t>(low_rank);
  const double low = samples[k];
  const double high = samples[k + 1];
  return low + (high - low) * dist;
}
// Winsorize the samples to clean up any very extreme outliers
// Will replace all samples below `percent` and above 100 - `percent` percentiles
// with the value at the percentile
// NOTE: Assumes the samples have been sorted, as we make use of PercentileValue
void Winsorize(std::vector<double> &samples, const double percent){
  const double low_percentile = PercentileValue(samples, percent);
  const double high_percentile = PercentileValue(samples, 100.0 - percent);
  for (std::vector<double>::iterator it = samples.begin(); it != samples.end(); ++it){
    if (*it < low_percentile){
      *it = low_percentile;
    }
    else if (*it > high_percentile){
      *it = high_percentile;
    }
  }
}
// Compute the mean value of the dataset
double Mean(const std::vector<double> &samples){
  double mean = 0;
  for (std::vector<double>::const_iterator it = samples.begin(); it != samples.end(); ++it){
    mean += *it;
  }
  return mean / static_cast<double>(samples.size());
}
// Compute the sample variance of the samples
double Variance(const std::vector<double> &samples){
  double mean = Mean(samples);
  double square_deviations = 0;
  for (std::vector<double>::const_iterator it = samples.begin(); it != samples.end(); ++it){
    square_deviations += std::pow(*it - mean, 2.0);
  }
  return square_deviations / (static_cast<double>(samples.size()) - 1.0);
}
// Compute the standard deviation of the samples
double StandardDeviation(const std::vector<double> &samples){
  return std::sqrt(Variance(samples));
}
// Compute the median absolute deviation of the dataset
double MedianAbsDeviation(const std::vector<double> &samples){
  std::vector<double> abs_deviations;
  abs_deviations.reserve(samples.size());
  const double median = PercentileValue(samples, 50.0);
  for (std::vector<double>::const_iterator it = samples.begin(); it != samples.end(); ++it){
    abs_deviations.push_back(std::abs(*it - median));
  }
  std::sort(abs_deviations.begin(), abs_deviations.end());
  return PercentileValue(abs_deviations, 50.0);
}
} // stats
