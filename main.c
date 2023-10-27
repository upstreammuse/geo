#include <assert.h>
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

double maxFinError = -1;
TestCaseDegree maxFinErrorCase;
double maxEndLatError = -1;
TestCaseDegree maxEndLatErrorCase;
double maxEndLonError = -1;
TestCaseDegree maxEndLonErrorCase;

void printEndpointResults() {
   printf("Worst final bearing error: %e\n", maxFinError);
   printTestCaseDegree(maxFinErrorCase);
   printf("\n");
   printf("Worst end latitude error: %e\n", maxEndLatError);
   printTestCaseDegree(maxEndLatErrorCase);
   printf("\n");
   printf("Worst end longitude error: %e\n", maxEndLonError);
   printTestCaseDegree(maxEndLonErrorCase);
   printf("\n");
}

void endpointTest(TestCaseDegree const tc) {
   TestCaseRadian workingR;
   tcToRadian(&workingR, tc);
   endpoint(&workingR.end, &workingR.rb, workingR.start);
   TestCaseDegree workingD;
   tcToDegree(&workingD, workingR);
   // best we can get convering to radians and back
   assert(fabs(tc.start.lat - workingD.start.lat) < 1e-13);
   assert(tc.start.lon == workingD.start.lon);
   assert(tc.rb.range == workingD.rb.range);
   // best we can get convering to radians and back
   assert(fabs(tc.rb.ini - workingD.rb.ini) < 1e-13);
   double error = fabs(tc.rb.fin - workingD.rb.fin);
   if (error > maxFinError) {
      maxFinError = error;
      maxFinErrorCase = tc;
   }
   error = fabs(tc.end.lat - workingD.end.lat);
   if (error > maxEndLatError) {
      maxEndLatError = error;
      maxEndLatErrorCase = tc;
   }
   error = fabs(tc.end.lon - workingD.end.lon);
   if (error > maxEndLonError) {
      maxEndLonError = error;
      maxEndLonErrorCase = tc;
   }
}

//const double thresholds[] = {
//   1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1, 1000, 1e100
//};
//int latcounts[sizeof (thresholds) / sizeof (double)] = {0};
//int loncounts[sizeof (thresholds) / sizeof (double)] = {0};
//int bearingcounts[sizeof (thresholds) / sizeof (double)] = {0};

int main(void) {
   char buffer[1000];
   FILE* f = fopen("GeodTest.dat", "rt");
//   double maxThreshold = -1;
//   TestCaseDegree maxTestCase;
   double maxThreshold_rb = -1;
   TestCaseDegree maxTestCase_rb;
   while (fgets(buffer, 999, f) != NULL) {
      /* read test case inputs */
      TestCaseDegree tc;
      sscanf(buffer, "%lf %lf %lf %lf %lf %lf %lf",
             &tc.start.lat, &tc.start.lon, &tc.rb.ini, &tc.end.lat, &tc.end.lon,
             &tc.rb.fin, &tc.rb.range);

//      /* calculate endpoint results */
//      TestCaseDegree myTC;
//      {
//         TestCaseRadian myTCR;
//         tcToRadian(&myTCR, tc);
//         endpoint(&myTCR.end, &myTCR.rb, myTCR.start);
//         tcToDegree(&myTC, myTCR);
//      }
      /* run endpoint test */
      endpointTest(tc);

      /* calculate range/bearing results */
      TestCaseDegree myTC_rb;
      {
         TestCaseRadian myTCR;
         tcToRadian(&myTCR, tc);
         rangeBearing(&myTCR.rb, myTCR.start, myTCR.end);
         tcToDegree(&myTC_rb, myTCR);
      }

//      /* save worst endpoint test case for comparison */
//      double threshold = fmax(fmax(fmax(fabs(myTC.end.lat - tc.end.lat),
//                                        fabs(myTC.end.lon - tc.end.lon)),
//                                   fabs(myTC.rb.ini - tc.rb.ini)),
//                              fabs(myTC.rb.fin - tc.rb.fin));
//      if (threshold > maxThreshold) {
//         maxTestCase = myTC;
//         maxThreshold = threshold;
//      }

      /* save worst range/bearing test case for comparison */
      double threshold = fmax(fmax(fabs(myTC_rb.rb.range - tc.rb.range),
                            fabs(myTC_rb.rb.ini - tc.rb.ini)),
                       fabs(myTC_rb.rb.fin - tc.rb.fin));
      if (threshold > maxThreshold_rb) {
         maxTestCase_rb = myTC_rb;
         maxThreshold_rb = threshold;
      }

//      /* increment the counter matching the threshold difference */
//      int i = 0;
//      while (fabs(myTC.end.lat - tc.end.lat) > thresholds[i]) i++;
//      latcounts[i]++;
//      i = 0;
//      while (fabs(myTC.end.lon - tc.end.lon) > thresholds[i]) i++;
//      loncounts[i]++;
//      i = 0;
//      while (fabs(myTC.rb.ini - tc.rb.ini) > thresholds[i]) i++;
//      bearingcounts[i]++;
   }
   fclose(f);

   printEndpointResults();

//   for (int i = 0; thresholds[i] != 1e100; i++) {
//      printf("%d lat, %d lon, %d bearing within %e\n", latcounts[i], 
//             loncounts[i], bearingcounts[i], thresholds[i]);
//   }
//   printf("Worst endpoint test case: ");
//   printTestCaseDegree(maxTestCase);
//   printf("\nOff by: %e\n", maxThreshold);
   printf("Worst range/bearing test case: ");
   printTestCaseDegree(maxTestCase_rb);
   printf("\nOff by: %e\n", maxThreshold_rb);
   return 0;
}
