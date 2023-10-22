#include <math.h>
#include <stdio.h>
#include "geo.h"

typedef struct {
   PositionDegree start;
   PositionDegree end;
   RangeBearingDegree rb;
} TestCaseDegree;

typedef struct {
   PositionRadian start;
   PositionRadian end;
   RangeBearingRadian rb;
} TestCaseRadian;

void tcToDegree(TestCaseDegree* out, TestCaseRadian in) {
   posToDegree(&out->start, in.start);
   posToDegree(&out->end, in.end);
   rbToDegree(&out->rb, in.rb);
}

void tcToRadian(TestCaseRadian* out, TestCaseDegree in) {
   posToRadian(&out->start, in.start);
   posToRadian(&out->end, in.end);
   rbToRadian(&out->rb, in.rb);
}

int printTestCaseDegree(TestCaseDegree tc) {
   return printf("(%lf %lf) via (%lf %lf) to (%lf %lf) arriving via %lf",
                 tc.start.lat, tc.start.lon, tc.rb.range, tc.rb.ini, tc.end.lat,
                 tc.end.lon, tc.rb.fin);
}

const double thresholds[] = {
   1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1, 1000, 1e100
};
int latcounts[sizeof (thresholds) / sizeof (double)] = {0};
int loncounts[sizeof (thresholds) / sizeof (double)] = {0};
int bearingcounts[sizeof (thresholds) / sizeof (double)] = {0};

int main(void) {
   char buffer[1000];
   FILE* f = fopen("GeodTest.dat", "rt");
   double maxThreshold = 0;
   TestCaseDegree maxTestCase;
   while (fgets(buffer, 999, f) != NULL) {
      /* read inputs */
      TestCaseDegree tc;
      sscanf(buffer, "%lf %lf %lf %lf %lf %lf %lf",
             &tc.start.lat, &tc.start.lon, &tc.rb.ini, &tc.end.lat, &tc.end.lon,
             &tc.rb.fin, &tc.rb.range);

      /* calculate results */
      TestCaseDegree myTC;
      {
         TestCaseRadian myTCR;
         tcToRadian(&myTCR, tc);
         endpoint(&myTCR.end, &myTCR.rb, myTCR.start);
         tcToDegree(&myTC, myTCR);
      }

      /* save worst test case for comparison */
      double threshold = fmax(fmax(fmax(fabs(myTC.end.lat - tc.end.lat),
                                        fabs(myTC.end.lon - tc.end.lon)),
                                   fabs(myTC.rb.ini - tc.rb.ini)),
                              fabs(myTC.rb.fin - tc.rb.fin));
      if (threshold > maxThreshold) {
         maxTestCase = myTC;
         maxThreshold = threshold;
      }

      /* increment the counter matching the threshold difference */
      int i = 0;
      while (fabs(myTC.end.lat - tc.end.lat) > thresholds[i]) i++;
      latcounts[i]++;
      i = 0;
      while (fabs(myTC.end.lon - tc.end.lon) > thresholds[i]) i++;
      loncounts[i]++;
      i = 0;
      while (fabs(myTC.rb.ini - tc.rb.ini) > thresholds[i]) i++;
      bearingcounts[i]++;
   }
   fclose(f);

   for (int i = 0; thresholds[i] != 1e100; i++) {
      printf("%d lat, %d lon, %d bearing within %e\n", latcounts[i], 
             loncounts[i], bearingcounts[i], thresholds[i]);
   }
   printf("Worst test case: ");
   printTestCaseDegree(maxTestCase);
   printf("\nOff by: %e\n", maxThreshold);
   return 0;
}
