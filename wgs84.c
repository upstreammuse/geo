#include "wgs84.h"

double const WGS84_A = 6378137.0;
double const WGS84_F = 1 / 298.257223563;

double WGS84_B() {
   return WGS84_A - WGS84_A * WGS84_F;
}
