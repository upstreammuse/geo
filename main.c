#include <math.h>
#include <stdio.h>
#include "geo.h"

int main(void) {
   char buffer[1000];
   FILE* f = fopen("GeodTest.dat", "rt");
   int fail = 0;
   int lines = 0;
   while (fgets(buffer, 999, f) != NULL) {
      lines++;
      double lat1, lon1, lat2, lon2, range, bearing1, bearing2;
      sscanf(buffer, "%lf %lf %lf %lf %lf %lf %lf",
             &lat1, &lon1, &bearing1, &lat2, &lon2, &bearing2, &range);
      lat1 = lat1 / 180 * M_PI;
      bearing1 = bearing1 / 180 * M_PI;
      double myLat2, myLon2, myBearing2;
      direct(&myLat2, &myLon2, &myBearing2, lat1, range, bearing1);
      myLat2 = myLat2 / M_PI * 180;
      myLon2 = myLon2 / M_PI * 180;
      /* TODO these are needed because abs results are > PI sometimes */
      if (myLon2 < -180) {
         myLon2 += 360;
      }
      if (myLon2 > 180) {
         myLon2 -= 360;
      }
      myBearing2 = myBearing2 / M_PI * 180;
      /* need to be correct to 3rd decimal digit of minutes */
      if (fabs(myLat2 - lat2) > 1.0 / (60 * 1000)) {
         printf("Theirs: %lf   Mine: %lf\n", lat2, myLat2);
         fail++;
      }
      if (fabs(myLon2 - lon2) > 1.0 / (60 * 1000)) {
         printf("Theirs: %lf   Mine: %lf\n", lon2, myLon2);
         fail++;
      }
   }
   printf("%d failures of %d tests\n", fail, lines);
   fclose(f);
   return 0;
}
