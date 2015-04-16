typedef struct vtkData
{
    double x,y,z,w;
}vtkDataArray;

typedef struct splineData {
    int     order;
    int     numPoints;
    vtkDataArray *controlPoints;
    float   *knots;
}SplineData;

SplineData *InsertKnot(SplineData *oldCurveData, float tNew)
{
   SplineData  *newCurveData;
   int k,n;
   vtkDataArray  *b;             // Old control point vector
   vtkDataArray  *bHat;          // New control point vector
   float *x;                    // Old knot vector
   float *xHat;          // New knot vector
   float alpha;          // Interpolation ratio
   unsigned long i;              // Knot to insert after
   unsigned long j;              // Knot index for search
   bool foundIndex;     // Insertion index found?

   // Set up local variables for readability.
   k = oldCurveData->order;
   n = oldCurveData->numPoints;
   x = oldCurveData->knots;
   b = oldCurveData->controlPoints;

   // Allocate space for new control points and knot vector.
   bHat = malloc((n + 1) * sizeof(vtkDataArray));
   xHat = malloc((n + k + 1) * sizeof(float));

   // Allocate data structure for new curve.
   newCurveData = malloc(sizeof(SplineData));
   newCurveData->order = k;
   newCurveData->numPoints = n + 1;
   newCurveData->controlPoints = bHat;
   newCurveData->knots = xHat;

   // Find where to insert the new knot.
   for (j = 0, foundIndex = false; j < n + k; j++) {
      if (tNew > x[j] && tNew <= x[j + 1]) {
         i = j;
         foundIndex = true;
         break;
      }
   }

   // Return if not found.
   if (!foundIndex) {
      return (NULL);
   }
   // Copy knots to new vector.
   for (j = 0; j < n + k + 1; j++) {
      if (j <= i) {
         xHat[j] = x[j];
      } else if (j == i + 1) {
         xHat[j] = tNew;
      } else {
         xHat[j] = x[j - 1];
      }
   }

   // Compute position of new control point and new positions of
   // existing ones.
   for (j = 0; j < n + 1; j++) {
      if (j <= i - k + 1) {
         alpha = 1;
      } else if (i - k + 2 <= j && j <= i) {
         if (x[j + k - 1] - x[j] == 0) {
            alpha = 0;
         } else {
            alpha = (tNew - x[j]) / (x[j + k - 1] - x[j]);
         }
      } else {
         alpha = 0;
      }

      if (alpha == 0) {
         bHat[j] = b[j - 1];
      } else if (alpha == 1) {
         bHat[j] = b[j];
      } else {
         bHat[j].x = (1 - alpha) * b[j - 1].x + alpha * b[j].x;
         bHat[j].y = (1 - alpha) * b[j - 1].y + alpha * b[j].y;
         bHat[j].z = (1 - alpha) * b[j - 1].z + alpha * b[j].z;
         bHat[j].w = (1 - alpha) * b[j - 1].w + alpha * b[j].w;
      }
   }

   return (newCurveData);
}
