#include <Rcpp.h>
#include "Contract.h"

Contract::Contract()
{
  strike = 100;
  maturity = 30.0 / 252.0;
  type = "call";
}
