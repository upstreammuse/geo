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

int inverse(double* s, double* alpha1, double* alpha2,
            double const phi1, double const phi2, double const L) {

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

      if (fabs(lambda) >= M_PI) {
         return -1;
//         printf("lambda is %e, switching to alternate method\n", lambda);
//         inverse2(s, alpha1, alpha2, phi1, phi2, L);
//         return;
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
   /* "The common term in 14 and 20 should be calculated once and reused." */
   *alpha1 = atan2(cosU2 * sinLambda, term14);

   /* 21 */
   *alpha2 = atan2(cosU1 * sinLambda, -sinU1 * cosU2 + cosU1 * sinU2 * cosLambda);

   return 0;
}

// from the followup paper
void inverse2(double* s, double* alpha1, double* alpha2,
              double const phi1, double const phi2, double const L) {

   /* sanity check inputs */
   assert(!isnan(phi1) && fabs(phi1) <= M_PI/2);
   assert(!isnan(phi2) && fabs(phi2) <= M_PI/2);
   assert(!isnan(L) && fabs(L) <= M_PI);

   if (phi1 == -phi2) {
      printf("Special case alert!\n");
   }

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
   double const* const U2set = reducedLat(phi2);
   double const sinU2 = U2set[0];
   double const cosU2 = U2set[1];
   double const tanU2 = U2set[2];
   double const U2 = asin(sinU2);

   // initial approximations and values
   const double Lprime = (L < 0) ? (-M_PI - L) : (M_PI - L);
   assert(!isnan(Lprime) && fabs(Lprime) <= M_PI);
   /* 6 */ double cos2SigmaM = 0;
   /* 8 */ double cosSqAlpha = 0.5;
   /* 8 */ double oldSinAlpha;
   /* 8 */ double sinAlpha = sqrt(1 - cosSqAlpha);
   /* 9 */ double lambdaPrime = 0;
   /* 10 */ double sigma = M_PI - fabs(U1 + U2);
   assert(!isnan(sigma) && fabs(sigma) <= M_PI);
   int iterCount = 0;
   do {
      /* 5 */
      printf("cosSqAlpha=%e\n", cosSqAlpha);
      assert(!isnan(cosSqAlpha) && fabs(sqrt(cosSqAlpha)) <= 1);
      const double C = f / 16 * cosSqAlpha * (4 + f * (4 - 3 * cosSqAlpha));
      printf("C=%e\n", C);
      assert(!isnan(C));

      /* 6 */
      printf("cosSqAlpha=%e\n", cosSqAlpha);
      assert(cosSqAlpha != 0);
      cos2SigmaM = cos(sigma) - 2 * sinU1 * sinU2 / cosSqAlpha;
      printf("cos2SM=%e\n", cos2SigmaM);
      if (cos2SigmaM < -1 || cos2SigmaM > 1) {
         printf("correcting cos2SigmaM out of range\n");
         if (cos2SigmaM < -1) cos2SigmaM = -1;
         if (cos2SigmaM > 1) cos2SigmaM = 1;
      }
      assert(!isnan(cos2SigmaM) && fabs(cos2SigmaM) <= 1);

      /* 7 */
      const double D = (1 - C) * f * (sigma + C * sin(sigma) * (cos2SigmaM + C * cos(sigma) * (-1 + 2 * cos2SigmaM * cos2SigmaM)));
      assert(!isnan(D));

      /* 8 */
      oldSinAlpha = sinAlpha;
      assert(!isnan(oldSinAlpha) && fabs(oldSinAlpha) <= 1);
      printf("L'=%e l'=%e D=%e\n", Lprime, lambdaPrime, D);
      sinAlpha = (Lprime - lambdaPrime) / D;
      printf("sinAlpha=%e\n", sinAlpha);
      assert(!isnan(sinAlpha) && fabs(sinAlpha) <= 1);
      cosSqAlpha = 1 - sinAlpha * sinAlpha;
      printf("cosSqAlpha=%e\n", cosSqAlpha);
      assert(!isnan(cosSqAlpha) && fabs(sqrt(cosSqAlpha)) <= 1);

      /* 9 */
      printf("sA=%e sS=%e cU1=%e cU2=%e\n", sinAlpha, sin(sigma), cosU1, cosU2);
      const double sinLambdaPrime = sinAlpha * sin(sigma) / (cosU1 * cosU2);
      printf("sinl'=%e\n", sinLambdaPrime);
      assert(!isnan(sinLambdaPrime) && fabs(sinLambdaPrime) <= 1);
      lambdaPrime = asin(sinLambdaPrime);
      printf("l'=%e\n", lambdaPrime);
      assert(!isnan(lambdaPrime));

      /* 10 */
      const double t1 = cosU2 * sinLambdaPrime;
      const double t2 = cosU1 * sinU2 + sinU1 * cosU2 * cos(lambdaPrime);
      const double sinSqSigma = t1 * t1 + t2 * t2;
      printf("sinSqSigma=%e\n", sinSqSigma);
      assert(!isnan(sinSqSigma) && fabs(sqrt(sinSqSigma)) <= 1);
      sigma = asin(sqrt(sinSqSigma));
      printf("sigma=%e\n", sigma);
      assert(!isnan(sigma) && fabs(sigma) <= M_PI);
   } while (fabs(sinAlpha - oldSinAlpha) > 1e-12 && iterCount++ < 1000);

   const double sinSigma = sin(sigma);
   const double cosSigma = cos(sigma);
   const double cosSq2SigmaM = cos2SigmaM * cos2SigmaM;
   const double sinSqSigma = sinSigma * sinSigma;

   /* 11 */
   const double sinAlpha1 = sinAlpha / cosU1;

   /* 12 */
   double cosAlpha1 = sqrt(1 - sinAlpha1 * sinAlpha1);

   /* 13 */
   if (cosU1 * sinU2 + sinU1 * cosU2 * cos(lambdaPrime) < 0) {
      cosAlpha1 = -cosAlpha1;
   }

   /* 14 */
   *alpha2 = atan2(sinAlpha, -sinU1 * sin(sigma) + cosU1 * cos(sigma) * cosAlpha1);
   assert(!isnan(*alpha2));

   /* additional commentary */
   *alpha1 = atan2(sinAlpha1, cosAlpha1);
   assert(!isnan(*alpha1));

   /* 15 */
   const double epsilon = (a * a - b * b) / (b * b);
//   printf("epsilon=%e\n", epsilon);
   const double E = sqrt(1 + epsilon * cosSqAlpha);
//   printf("E=%e\n", E);

   /* 16 */
   const double F = (E - 1) / (E + 1);
//   printf("F=%e\n", F);

   /* 17 */
   const double A = (1 + F * F / 4) / (1 - F);
//   printf("A=%e\n", A);

   /* 18 */
   const double B = F * (1 - F * F * 3 / 8);

   /* 19 */
   const double deltaSigma = B * sinSigma * (cos2SigmaM + B / 4 * (cosSigma * (-1 + 2 * cosSq2SigmaM) - B / 6 * cos2SigmaM * (-3 + 4 * sinSqSigma) * (-3 + 4 * cosSq2SigmaM)));

   /* 20 */
   *s = b * A * (sigma - deltaSigma);
   assert(!isnan(*s));
//   printf("%e\n", *s);
}
