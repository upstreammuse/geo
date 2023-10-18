#include "geo.h"

#include <math.h>

void direct(double* phi2, double* L, double* alpha2, double phi1, double s, double alpha1) {
   /* WGS-84 definitions */
   const double a = 6378137.0;
   const double f = 1 / 298.257223563;
   const double b = a - a * f;

   const double tanU1 = (1 - f) * tan(phi1);
   /* OK to do this because ph1 is 4th and 1st quadrant, so don't need negative
    cosine values */
   const double cosU1 = 1 / sqrt(1 + tanU1 * tanU1);
   const double sinU1 = tanU1 * cosU1;

   /* (1) */
   const double cosAlpha1 = cos(alpha1);
   const double sigma1 = atan2(tanU1, cosAlpha1);

   /* (2) */
   const double sinAlpha1 = sin(alpha1);
   const double sinAlpha = cosU1 * sinAlpha1;
   const double sinSqAlpha = sinAlpha * sinAlpha;

   /* 3 */
   const double cosSqAlpha = 1 - sinSqAlpha;
   const double uSq = cosSqAlpha * (a * a - b * b) / (b * b);
   const double A = 1 + uSq / 16384 * (4096 + uSq * (-768 + uSq * (320 - 175 * uSq)));

   /* 4 */
   const double B = uSq / 1024 * (256 + uSq * (-128 + uSq * (74 - 47 * uSq)));

   /* use first term of 7 to initialize sigma */
   double sigma = s / (b * A);
   double oldSigma = 0;

   double sinSigma;
   double cosSigma;
   double cosTwoSigmaM;
   double cosSqTwoSigmaM;
   do {
      /* 5 */
      const double twoSigmaM = 2 * sigma1 + sigma;

      /* 6 */
      sinSigma = sin(sigma);
      cosSigma = cos(sigma);
      const double sinSqSigma = sinSigma * sinSigma;
      cosTwoSigmaM = cos(twoSigmaM);
      cosSqTwoSigmaM = cosTwoSigmaM * cosTwoSigmaM;
      const double deltaSigma = B * sinSigma * (cosTwoSigmaM + .25 * B * (cosSigma * (-1 + 2 * cosSqTwoSigmaM) - 1.0/6.0 * B * cosTwoSigmaM * (-3 + 4 * sinSqSigma) * (-3 + 4 * cosSqTwoSigmaM)));

      /* 7 */
      oldSigma = sigma;
      sigma = s / (b * A) + deltaSigma;

      /* Vincenty's suggested threshold */
   } while (fabs(sigma - oldSigma) > 1e-12);

   /* 8 */
   double term8and12 = sinU1 * sinSigma - cosU1 * cosSigma * cosAlpha1;
   *phi2 = atan2(sinU1 * cosSigma + cosU1 * sinSigma * cosAlpha1,
                 (1 - f) * sqrt(sinSqAlpha + term8and12 * term8and12));

   /* 9 */
   const double lambda = atan2(sinSigma * sinAlpha1,
                               cosU1 * cosSigma - sinU1 * sinSigma * cosAlpha1);

   /* 10 */
   const double C = f / 16 * cosSqAlpha * (4 + f * (4 - 3 * cosSqAlpha));

   /* 11 */
   *L = lambda - (1 - C) * f * sinAlpha * (sigma + C * sinSigma * (cosTwoSigmaM + C * cosSigma * (-1 + 2 * cosSqTwoSigmaM)));

   /* 12 */
   *alpha2 = atan2(sinAlpha, -term8and12);
}
