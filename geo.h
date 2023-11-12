#pragma once

typedef struct {
   double lat;
   double lon;
} PositionDegree;

typedef struct {
   double lat;
   double lon;
} PositionRadian;

typedef struct {
   double range;
   double ini;
   double fin;
} RangeBearingDegree;

typedef struct {
   double range;
   double ini;
   double fin;
} RangeBearingRadian;

void posToDegree(PositionDegree*, PositionRadian);
void posToRadian(PositionRadian*, PositionDegree);
void rbToDegree(RangeBearingDegree*, RangeBearingRadian);
void rbToRadian(RangeBearingRadian*, RangeBearingDegree);

void endpoint(PositionRadian* end, RangeBearingRadian*, PositionRadian start);
int rangeBearing(RangeBearingRadian*, PositionRadian start, PositionRadian end);

void direct(double* phi2, double* L, double* alpha2, double phi1, double s, double alpha1);
int inverse(double* s, double* alpha1, double* alpha2, double phi1, double phi2, double L);
int inverse2(double* s, double* alpha1, double* alpha2, double phi1, double phi2, double L);
