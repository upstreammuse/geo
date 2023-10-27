#include "geo.h"

#include <math.h>

#define DEGREES(rad) (rad / M_PI * 180)
#define RADIANS(deg) (deg / 180 * M_PI)

void posToDegree(PositionDegree* out, PositionRadian in) {
   out->lat = DEGREES(in.lat);
   out->lon = DEGREES(in.lon);
}

void posToRadian(PositionRadian* out, PositionDegree in) {
   out->lat = RADIANS(in.lat);
   out->lon = RADIANS(in.lon);
}

void rbToDegree(RangeBearingDegree* out, RangeBearingRadian in) {
   out->range = in.range;
   out->ini = DEGREES(in.ini);
   out->fin = DEGREES(in.fin);
}

void rbToRadian(RangeBearingRadian* out, RangeBearingDegree in) {
   out->range = in.range;
   out->ini = RADIANS(in.ini);
   out->fin = RADIANS(in.fin);
}

void endpoint(PositionRadian* end, RangeBearingRadian* rb,
              PositionRadian start) {
   direct(&end->lat, &end->lon, &rb->fin, start.lat, rb->range, rb->ini);
   end->lon += start.lon;
   if (end->lon < -M_PI) {
      end->lon += 2 * M_PI;
   }
   if (end->lon > M_PI) {
      end->lon -= 2 * M_PI;
   }
   if (rb->fin < 0) {
      rb->fin += 2 * M_PI;
   }
}

void rangeBearing(RangeBearingRadian* rb, PositionRadian start, PositionRadian end) {
   double L = end.lon - start.lon;
   if (L < -M_PI) {
      L += 2 * M_PI;
   }
   if (L > M_PI) {
      L -= 2 * M_PI;
   }
   inverse2(&rb->range, &rb->ini, &rb->fin, start.lat, end.lat, L);
}

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
      const double deltaSigma = B * sinSigma * (cosTwoSigmaM + B / 4 * (cosSigma * (-1 + 2 * cosSqTwoSigmaM) - B / 6 * cosTwoSigmaM * (-3 + 4 * sinSqSigma) * (-3 + 4 * cosSqTwoSigmaM)));

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

void inverse(double* s, double* alpha1, double* alpha2, double phi1, double phi2, double L) {
   /* WGS-84 definitions */
   const double a = 6378137.0;
   const double f = 1 / 298.257223563;
   const double b = a - a * f;

   /* setup */
   const double tanU1 = (1 - f) * tan(phi1);
   const double tanU2 = (1 - f) * tan(phi2);
   /* OK to do this because ph1 is 4th and 1st quadrant, so don't need negative
    cosine values */
   const double cosU1 = 1 / sqrt(1 + tanU1 * tanU1);
   const double cosU2 = 1 / sqrt(1 + tanU2 * tanU2);
   const double sinU1 = tanU1 * cosU1;
   const double sinU2 = tanU2 * cosU2;

   /* 13 */
   double lambda = L;

   double sinLambda;
   double cosLambda;
   double term14;
   double oldLambda;
   double sigma;
   double cosSq2SigmaM;
   double sinSqSigma;
   double cosSigma;
   double cos2SigmaM;
   double sinSigma;
   double cosSqAlpha;
   double iterCount = 0;
   do {
      /* 14 */
      sinLambda = sin(lambda);
      cosLambda = cos(lambda);
      term14 = cosU1 * sinU2 - sinU1 * cosU2 * cosLambda;
      sinSqSigma = (cosU2 * sinLambda) * (cosU2 * sinLambda) + term14 * term14;

      /* 15 */
      cosSigma = sinU1 * sinU2 + cosU1 * cosU2 * cosLambda;

      /* 16 */
      // TODO is the loss of negative values an issue here?
      sinSigma = sqrt(sinSqSigma);
      sigma = atan2(sinSigma, cosSigma);

      /* 17 */
      double sinAlpha = cosU1 * cosU2 * sinLambda / sinSigma;

      /* 18 */
      cosSqAlpha = 1 - sinAlpha * sinAlpha;
      cos2SigmaM = cosSigma - 2 * sinU1 * sinU2 / cosSqAlpha;
      cosSq2SigmaM = cos2SigmaM * cos2SigmaM;

      /* 10 */
      const double C = f / 16 * cosSqAlpha * (4 + f * (4 - 3 * cosSqAlpha));

      /* 11 */
      oldLambda = lambda;
      lambda = L + (1 - C) * f * sinAlpha * (sigma + C * sinSigma * (cos2SigmaM + C * cosSigma * (-1 + 2 * cosSq2SigmaM)));
   } while (fabs(lambda - oldLambda) > 1e-12 && iterCount++ < 1000);

   /* 3 */
   const double uSq = cosSqAlpha * (a * a - b * b) / (b * b);
   const double A = 1 + uSq / 16384 * (4096 + uSq * (-768 + uSq * (320 - 175 * uSq)));

   /* 4 */
   const double B = uSq / 1024 * (256 + uSq * (-128 + uSq * (74 - 47 * uSq)));

   /* 6 */
   const double deltaSigma = B * sinSigma * (cos2SigmaM + B / 4 * (cosSigma * (-1 + 2 * cosSq2SigmaM) - B / 6 * cos2SigmaM * (-3 + 4 * sinSqSigma) * (-3 + 4 * cosSq2SigmaM)));

   /* 19 */
   *s = b * A * (sigma - deltaSigma);

   /* 20 */
   /* "The common term in 14 and 20 should be calculated once and reused." */
   *alpha1 = atan2(cosU2 * sinLambda, term14);

   /* 21 */
   *alpha2 = atan2(cosU1 * sinLambda, -sinU1 * cosU2 + cosU1 * sinU2 * cosLambda);
}

// from the followup paper
void inverse2(double* s, double* alpha1, double* alpha2, double phi1, double phi2, double L) {
   /* WGS-84 definitions */
   const double a = 6378137.0;
   const double f = 1 / 298.257223563;
   const double b = a - a * f;

   // definitions
   const double tanU1 = (1 - f) * tan(phi1);
   const double tanU2 = (1 - f) * tan(phi2);
   /* OK to do this because ph1 is 4th and 1st quadrant, so don't need negative
    cosine values */
   const double cosU1 = 1 / sqrt(1 + tanU1 * tanU1);
   const double cosU2 = 1 / sqrt(1 + tanU2 * tanU2);
   const double sinU1 = tanU1 * cosU1;
   const double sinU2 = tanU2 * cosU2;
   const double U1 = asin(sinU1);
   const double U2 = asin(sinU2);

   // initial approximations and values
   const double Lprime = L < 0 ? -M_PI - L : M_PI - L;
   double lambdaPrime = 0;
   double cosSqAlpha = 0.5;
   double cos2SigmaM = 0;
   double sigma = M_PI - fabs(U1 + U2);

   double sinAlpha = sqrt(0.5);  // based on initial definition of cosSqAlpha
   double oldSinAlpha;
   double cosLambdaPrime;
   double sinSigma;
   double cosSigma;
   double cosSq2SigmaM;
   double sinSqSigma;
   do {
      /* 5 */
      const double C = f / 16 * cosSqAlpha * (4 + f * (4 - 3 * cosSqAlpha));

      /* 6 */
      cosSigma = cos(sigma);
      cos2SigmaM = cosSigma - 2 * sinU1 * sinU2 / cosSqAlpha;

      /* 7 */
      cosSq2SigmaM = cos2SigmaM * cos2SigmaM;
      sinSigma = sin(sigma);
      sigma = asin(sigma);
      const double D = (1 - C) * f * (sigma + C * sinSigma * (cos2SigmaM + C * cosSigma * (-1 + 2 * cosSq2SigmaM)));

      /* 8 */
      oldSinAlpha = sinAlpha;
      sinAlpha = (Lprime - lambdaPrime) / D;
      cosSqAlpha = 1 - sinAlpha * sinAlpha;

      /* 9 */
      const double sinLambdaPrime = sinAlpha * sinSigma / (cosU1 * cosU2);
      lambdaPrime = asin(sinLambdaPrime);

      /* 10 */
      // TODO loss of sign info
      cosLambdaPrime = sqrt(1 - sinLambdaPrime * sinLambdaPrime);
      const double t1 = cosU2 * sinLambdaPrime;
      const double t2 = cosU1 * sinU2 + sinU1 * cosU2 * cosLambdaPrime;
      sinSqSigma = t1 * t1 + t2 * t2;
   } while (fabs(sinAlpha - oldSinAlpha) > 1e-12);

   /* 11 */
   const double sinAlpha1 = sinAlpha / cosU1;

   /* 12 */
   double cosAlpha1 = sqrt(1 - sinAlpha1 * sinAlpha1);

   /* 13 */
   if (cosU1 * sinU2 + sinU1 * cosU2 * cosLambdaPrime < 0) {
      cosAlpha1 = -cosAlpha1;
   }

   /* 14 */
   *alpha2 = atan2(sinAlpha, -sinU1 * sinSigma + cosU1 * cosSigma * cosAlpha1);

   /* additional commentary */
   *alpha1 = atan2(sinAlpha1, cosAlpha1);

   /* 15 */
   const double epsilon = (a * a - b * b) / (b * b);
   const double E = sqrt(1 + epsilon * cosSqAlpha);

   /* 16 */
   const double F = (E - 1) / (E + 1);

   /* 17 */
   const double A = (1 + F * F / 4) / (1 - F);

   /* 18 */
   const double B = F * (1 - F * F * 3 / 8);

   /* 19 */
   const double deltaSigma = B * sinSigma * (cos2SigmaM + B / 4 * (cosSigma * (-1 + 2 * cosSq2SigmaM) - B / 6 * cos2SigmaM * (-3 + 4 * sinSqSigma) * (-3 + 4 * cosSq2SigmaM)));

   /* 20 */
   *s = b * A * (sigma - deltaSigma);
}
