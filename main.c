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

struct {
   double maxFinError;
   TestCaseDegree maxFinErrorCase;
   double maxEndLatError;
   TestCaseDegree maxEndLatErrorCase;
   double maxEndLonError;
   TestCaseDegree maxEndLonErrorCase;
} endpointResults = {.maxFinError = -1, .maxEndLatError = -1, .maxEndLonError = -1};

void printEndpointResults() {
   printf("===Endpoint Results===\n");
   printf("Worst final bearing error: %e\n", endpointResults.maxFinError);
   printTestCaseDegree(endpointResults.maxFinErrorCase);
   printf("\n");
   printf("Worst end latitude error: %e\n", endpointResults.maxEndLatError);
   printTestCaseDegree(endpointResults.maxEndLatErrorCase);
   printf("\n");
   printf("Worst end longitude error: %e\n", endpointResults.maxEndLonError);
   printTestCaseDegree(endpointResults.maxEndLonErrorCase);
   printf("\n");
}

void endpointTest(TestCaseDegree const tc) {
   /* run the test case and convert as needed */
   TestCaseRadian workingR;
   tcToRadian(&workingR, tc);
   endpoint(&workingR.end, &workingR.rb, workingR.start);
   TestCaseDegree workingD;
   tcToDegree(&workingD, workingR);

   /* check that the original values are still good */
   // best we can get convering to radians and back?
   assert(fabs(tc.start.lat - workingD.start.lat) < 1e-13);
   assert(tc.start.lon == workingD.start.lon);
   assert(tc.rb.range == workingD.rb.range);
   // best we can get convering to radians and back?
   assert(fabs(tc.rb.ini - workingD.rb.ini) < 1e-13);

   /* compare results to test cases */
   double error = fabs(tc.rb.fin - workingD.rb.fin);
   if (error > endpointResults.maxFinError) {
      endpointResults.maxFinError = error;
      endpointResults.maxFinErrorCase = tc;
   }
   error = fabs(tc.end.lat - workingD.end.lat);
   if (error > endpointResults.maxEndLatError) {
      endpointResults.maxEndLatError = error;
      endpointResults.maxEndLatErrorCase = tc;
   }
   error = fabs(tc.end.lon - workingD.end.lon);
   if (error > endpointResults.maxEndLonError) {
      endpointResults.maxEndLonError = error;
      endpointResults.maxEndLonErrorCase = tc;
   }
}

struct {
   double maxIniError;
   TestCaseDegree maxIniErrorCase;
   double maxFinError;
   TestCaseDegree maxFinErrorCase;
   double maxRangeError;
   TestCaseDegree maxRangeErrorCase;
   int numfails;
   int numtests;
} rangeBearingResults = {.maxIniError = -1, .maxFinError = -1, .maxRangeError = -1, .numfails = 0, .numtests = 0};

void printRangeBearingResults() {
   printf("===Range Bearing Results===\n");
   printf("Worst initial bearing error: %e\n", rangeBearingResults.maxIniError);
   printTestCaseDegree(rangeBearingResults.maxIniErrorCase);
   printf("\n");
   printf("Worst final bearing error: %e\n", rangeBearingResults.maxFinError);
   printTestCaseDegree(rangeBearingResults.maxFinErrorCase);
   printf("\n");
   printf("Worst range error: %e\n", rangeBearingResults.maxRangeError);
   printTestCaseDegree(rangeBearingResults.maxRangeErrorCase);
   printf("\n");
   printf("%d failed test cases out of %d\n", rangeBearingResults.numfails, rangeBearingResults.numtests);
}

void rangeBearingTest(TestCaseDegree const tc) {
   /* run the test case and convert as needed */
   TestCaseRadian workingR;
   tcToRadian(&workingR, tc);
   rangeBearingResults.numtests++;
   if (rangeBearing(&workingR.rb, workingR.start, workingR.end)) {
      rangeBearingResults.numfails++;
   }
   TestCaseDegree workingD;
   tcToDegree(&workingD, workingR);

   /* make sure original values are still good */
   // best we can get convering to radians and back?
   assert(fabs(tc.start.lat - workingD.start.lat) < 1e-13);
   assert(tc.start.lon == workingD.start.lon);
   // best we can get convering to radians and back?
   assert(fabs(tc.end.lat - workingD.end.lat) < 1e-13);
   assert(fabs(tc.end.lon - workingD.end.lon) < 1e-13);

   /* compare results to test case */
   double error = fabs(tc.rb.ini - workingD.rb.ini);
   if (error > rangeBearingResults.maxIniError) {
      rangeBearingResults.maxIniError = error;
      rangeBearingResults.maxIniErrorCase = tc;
   }
   error = fabs(tc.rb.fin - workingD.rb.fin);
   if (error > rangeBearingResults.maxFinError) {
      rangeBearingResults.maxFinError = error;
      rangeBearingResults.maxFinErrorCase = tc;
   }
   error = fabs(tc.rb.range - workingD.rb.range);
   if (error > rangeBearingResults.maxRangeError) {
      rangeBearingResults.maxRangeError = error;
      rangeBearingResults.maxRangeErrorCase = tc;
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
   while (fgets(buffer, 999, f) != NULL) {
      TestCaseDegree tc;
      sscanf(buffer, "%lf %lf %lf %lf %lf %lf %lf",
             &tc.start.lat, &tc.start.lon, &tc.rb.ini, &tc.end.lat, &tc.end.lon,
             &tc.rb.fin, &tc.rb.range);
      endpointTest(tc);
      rangeBearingTest(tc);
   }
   fclose(f);
   printEndpointResults();
   printRangeBearingResults();
   return 0;
}
