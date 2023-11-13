#include "geo.h"

#include <math.h>
#include <stdio.h>
#include <assert.h>

#include "wgs84.h"

#define DEGREES(rad) (rad / M_PI * 180)
#define RADIANS(deg) (deg / 180 * M_PI)

void posToDegree(PositionDegree* out, PositionRadian const in) {
   out->lat = DEGREES(in.lat);
   out->lon = DEGREES(in.lon);
}

void posToRadian(PositionRadian* out, PositionDegree const in) {
   out->lat = RADIANS(in.lat);
   out->lon = RADIANS(in.lon);
}

void rbToDegree(RangeBearingDegree* out, RangeBearingRadian const in) {
   out->range = in.range;
   out->ini = DEGREES(in.ini);
   out->fin = DEGREES(in.fin);
}

void rbToRadian(RangeBearingRadian* out, RangeBearingDegree const in) {
   out->range = in.range;
   out->ini = RADIANS(in.ini);
   out->fin = RADIANS(in.fin);
}

void endpoint(PositionRadian* end, RangeBearingRadian* rb,
              PositionRadian const start) {
   direct(&end->lat, &end->lon, &rb->fin, start.lat, rb->range, rb->ini);
   end->lon += start.lon;
   while (end->lon < -M_PI) {
      end->lon += 2 * M_PI;
   }
   while (end->lon > M_PI) {
      end->lon -= 2 * M_PI;
   }
   while (rb->fin < 0) {
      rb->fin += 2 * M_PI;
   }
}

int rangeBearing(RangeBearingRadian* rb, 
                 PositionRadian const start, PositionRadian const end) {
   double L = end.lon - start.lon;
//   if (L < -M_PI) {
//      L += 2 * M_PI;
//   }
//   if (L > M_PI) {
//      L -= 2 * M_PI;
//   }
   return inverse(&rb->range, &rb->ini, &rb->fin, start.lat, end.lat, L);
}

double* reducedLat(double const phi) {
   static double retval[3];
   double const f = WGS84_F;

   /* tanU */
   retval[2] = (1 - f) * tan(phi);
   assert(!isnan(retval[2]));

   /* cosU */
   /* OK to do this because ph1 is 4th and 1st quadrant, so don't need negative
    cosine values */
   retval[1] = 1 / sqrt(1 + retval[2] * retval[2]);
   assert(!isnan(retval[1]) && retval[1] >= 0 && retval[1] <= 1);

   /* sinU */
   retval[0] = retval[1] * retval[2];
   assert(!isnan(retval[0]) && fabs(retval[0]) <= 1);

   return retval;
}

void direct(double* phi2, double* L, double* alpha2, 
            double const phi1, double const s, double const alpha1) {

   /* WGS-84 definitions */
   const double a = WGS84_A;
   const double f = WGS84_F;
   const double b = WGS84_B();

   /* reduced latitudes */
   double const* const U1 = reducedLat(phi1);
   double const sinU1 = U1[0];
   double const cosU1 = U1[1];
   double const tanU1 = U1[2];

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
   const double term8and12 = sinU1 * sinSigma - cosU1 * cosSigma * cosAlpha1;
   *phi2 = atan2(sinU1 * cosSigma + cosU1 * sinSigma * cosAlpha1,
                 (1 - f) * sqrt(sinSqAlpha + term8and12 * term8and12));
   assert(!isnan(*phi2) && fabs(*phi2) <= M_PI / 2);

   /* 9 */
   const double lambda = atan2(sinSigma * sinAlpha1,
                               cosU1 * cosSigma - sinU1 * sinSigma * cosAlpha1);

   /* 10 */
   const double C = f / 16 * cosSqAlpha * (4 + f * (4 - 3 * cosSqAlpha));

   /* 11 */
   *L = lambda - (1 - C) * f * sinAlpha * (sigma + C * sinSigma * (cosTwoSigmaM + C * cosSigma * (-1 + 2 * cosSqTwoSigmaM)));
   assert(!isnan(*L));

   /* 12 */
   *alpha2 = atan2(sinAlpha, -term8and12);
   assert(!isnan(*alpha2) && *alpha2 >= 0 && *alpha2 <= 2 * M_PI);
}

int inverse(double* s, double* alpha1, double* alpha2,
            double const phi1, double const phi2, double const L) {

   /* check inputs */
   assert(!isnan(phi1) && fabs(phi1) <= M_PI / 2);
   assert(!isnan(phi2) && fabs(phi2) <= M_PI / 2);
   assert(!isnan(L) && fabs(L) <= M_PI);

   /* WGS-84 definitions */
   const double a = WGS84_A;
   const double f = WGS84_F;
   const double b = WGS84_B();

   /* reduced latitudes */
   double const* const U1 = reducedLat(phi1);
   double const sinU1 = U1[0];
   double const cosU1 = U1[1];
   double const tanU1 = U1[2];
   double const* const U2 = reducedLat(phi2);
   double const sinU2 = U2[0];
   double const cosU2 = U2[1];
   double const tanU2 = U2[2];

   /* 13 */
   double lambda = L;

   /* common */ double sinLambda;
   /* common */ double cosLambda;
   /* 14 */ double term14and20;
   /* 14 */ double sinSqSigma;
   /* 14 */ double sinSigma;
   /* 15 */ double cosSigma;
   /* 16 */ double sigma;
   /* 17 */ double cosSqAlpha;
   /* 18 */ double cos2SigmaM;
   /* 18 */ double cosSq2SigmaM;
   /* 11 */ double oldLambda;
   double iterCount = 0;
   do {
      /* common */
      sinLambda = sin(lambda);
      cosLambda = cos(lambda);

      /* 14 */
      term14and20 = cosU1 * sinU2 - sinU1 * cosU2 * cosLambda;
      sinSqSigma = (cosU2 * sinLambda) * (cosU2 * sinLambda) + term14and20 * term14and20;
      sinSigma = sqrt(sinSqSigma);

      /* 15 */
      cosSigma = sinU1 * sinU2 + cosU1 * cosU2 * cosLambda;

      /* 16 */
      sigma = atan2(sinSigma, cosSigma);

      /* 17 */
      double const sinAlpha = cosU1 * cosU2 * sinLambda / sinSigma;
      cosSqAlpha = 1 - sinAlpha * sinAlpha;

      /* 18 */
      cos2SigmaM = cosSigma - 2 * sinU1 * sinU2 / cosSqAlpha;
      cosSq2SigmaM = cos2SigmaM * cos2SigmaM;

      /* 10 */
      double const C = f / 16 * cosSqAlpha * (4 + f * (4 - 3 * cosSqAlpha));

      /* 11 */
      oldLambda = lambda;
      lambda = L + (1 - C) * f * sinAlpha * (sigma + C * sinSigma * (cos2SigmaM + C * cosSigma * (-1 + 2 * cosSq2SigmaM)));

      if (fabs(lambda) >= M_PI) {
         return inverse2(s, alpha1, alpha2, phi1, phi2, L);
      }

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
   *alpha1 = atan2(cosU2 * sinLambda, term14and20);

   /* 21 */
   *alpha2 = atan2(cosU1 * sinLambda, -sinU1 * cosU2 + cosU1 * sinU2 * cosLambda);

   return 0;
}

int sign(double x) {
   if (x < 0) return -1;
   if (x > 0) return 1;
   return 0;
}

// from the followup paper
int inverse2(double* s, double* alpha1, double* alpha2,
             double const phi1, double const phi2, double const L) {

   /* sanity check inputs */
   assert(!isnan(phi1) && fabs(phi1) <= M_PI/2);
   assert(!isnan(phi2) && fabs(phi2) <= M_PI/2);
   assert(!isnan(L) && fabs(L) <= M_PI);

   /* WGS-84 definitions */
   const double a = WGS84_A;
   const double f = WGS84_F;
   const double b = WGS84_B();

   /* reduced latitudes */
   double const* const U1set = reducedLat(phi1);
   double const sinU1 = U1set[0];
   double const cosU1 = U1set[1];
   double const tanU1 = U1set[2];
   double const U1 = asin(sinU1);
   assert(sign(U1) == sign(phi1));
   double const* const U2set = reducedLat(phi2);
   double const sinU2 = U2set[0];
   double const cosU2 = U2set[1];
   double const tanU2 = U2set[2];
   double const U2 = asin(sinU2);
   assert(sign(U2) == sign(phi2));

   // initial approximations and values
   const double Lprime = (L > 0) ? (M_PI - L) : (-M_PI - L);
   assert(!isnan(Lprime) && fabs(Lprime) <= M_PI);
   /* 6 */ double cos2SigmaM = 0;
   /* 8 */ double cosSqAlpha = 0.5;
   /* 8 */ double oldSinAlpha;
   /* 8 */ double sinAlpha = 1000;  /* obviously out of range, for first loop */
   /* 9 */ double lambdaPrime = 0;
   /* 10 */ double sigma = M_PI - fabs(U1 + U2);
   assert(!isnan(sigma) && sigma >= M_PI / 2 && sigma <= M_PI);
   int iterCount = 0;
   do {
      /* 5 */
      const double C = f / 16 * cosSqAlpha * (4 + f * (4 - 3 * cosSqAlpha));

      /* 6 */
      assert(cosSqAlpha != 0);
      cos2SigmaM = cos(sigma) - 2 * sinU1 * sinU2 / cosSqAlpha;
      /* per the paper, can't assume anything about this value */
      assert(!isnan(cos2SigmaM));

      /* 7 */
      const double D = (1 - C) * f * (sigma + C * sin(sigma) * (cos2SigmaM + C * cos(sigma) * (-1 + 2 * cos2SigmaM * cos2SigmaM)));
      assert(!isnan(D));

      /* 8 */
      oldSinAlpha = sinAlpha;
      sinAlpha = (Lprime - lambdaPrime) / D;
      cosSqAlpha = 1 - sinAlpha * sinAlpha;
      if (isnan(cosSqAlpha) || sqrt(cosSqAlpha) > 1) return -1;

      /* 9 */
      const double sinLambdaPrime = sinAlpha * sin(sigma) / (cosU1 * cosU2);
      if (isnan(sinLambdaPrime) || fabs(sinLambdaPrime) > 1) return -1;
      lambdaPrime = asin(sinLambdaPrime);
      assert(!isnan(lambdaPrime));

      /* 10 */
      const double t1 = cosU2 * sinLambdaPrime;
      const double t2 = cosU1 * sinU2 + sinU1 * cosU2 * cos(lambdaPrime);
      const double sinSqSigma = t1 * t1 + t2 * t2;
      assert(!isnan(sinSqSigma) && fabs(sqrt(sinSqSigma)) <= 1);
      sigma = asin(sqrt(sinSqSigma));
      assert(!isnan(sigma) && sigma >= 0 && sigma <= M_PI / 2);
      /* sigma is an antipodal angle in Q2, and needs mirroring from Q1 results */
      sigma = M_PI - sigma;
      assert(!isnan(sigma) && sigma >= M_PI / 2 && sigma <= M_PI);
   } while (fabs(sinAlpha - oldSinAlpha) > 1e-12 && iterCount++ < 1000);

   const double sinSigma = sin(sigma);
   const double cosSigma = cos(sigma);
   const double cosSq2SigmaM = cos2SigmaM * cos2SigmaM;
   const double sinSqSigma = sinSigma * sinSigma;

   /* 11 */
   // FIXME? this should not go out of range like this
   const double sinAlpha1 = (sinAlpha / cosU1) > 1 ? 1 : (sinAlpha / cosU1);

   /* 12 */
   double cosAlpha1 = sqrt(1 - sinAlpha1 * sinAlpha1);

   /* 13 */
   if (cosU1 * sinU2 + sinU1 * cosU2 * cos(lambdaPrime) < 0) {
      cosAlpha1 = -cosAlpha1;
   }

   /* 14 */
   *alpha2 = atan2(sinAlpha, -sinU1 * sin(sigma) + cosU1 * cos(sigma) * cosAlpha1);
   assert(!isnan(*alpha2) && *alpha1 >= 0 && *alpha1 <= 2 * M_PI);

   /* additional commentary */
   *alpha1 = atan2(sinAlpha1, cosAlpha1);
   assert(!isnan(*alpha1) && *alpha1 >= 0 && *alpha1 <= 2 * M_PI);

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
   assert(!isnan(*s));

   return 0;
}
