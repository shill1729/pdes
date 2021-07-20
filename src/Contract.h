#ifndef CONTRACT_H
#define CONTRACT_H

#include <Rcpp.h>
class Contract
{
public:
  Contract();
  Contract(double K, double T, std::string type);
  double strike;
  double maturity;
  std::string type;
  double payoff(double s);
private:
  double putPayoff(double s);
  double callPayoff(double s);
};
#endif
