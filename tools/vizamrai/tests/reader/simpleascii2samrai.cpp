//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/tools/vizamrai/tests/reader/simpleascii2samrai.cpp $
// Package:     Vizamrai
// Copyright:   (c) 1997-2000 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: A utility to read a Vizamrai output file and dump values
//              to the screen
//

#include <stdio.h>
#include <iostream.h>

#include <stdlib.h>

#include <rpc/types.h>
#include <rpc/xdr.h>

#define MAX_STRING 2056
#define SAMRAI_NULL_VALUE 1.0E30

#define VTK_LARGE_FLOAT 1.0e+38F
#define VTK_FLOAT_MIN -VTK_LARGE_FLOAT
#define VTK_FLOAT_MAX VTK_LARGE_FLOAT

int xdr_read_string(XDR *stream, char **data, int max_length)
{
   if (!xdr_string(stream, data, max_length))
   {
      return 1;
   }
   return 0;
}

// xdr read routines from Scott 
int xdr_read_int(XDR *stream, int *data, const int n)
{
   if (!xdr_vector(stream, (char *) data, n,
                   sizeof(int), (xdrproc_t) xdr_int))
   {
      return 1;
   }
   return 0;
}

int xdr_read_float(XDR *stream, float *data, const int n)
{
   if (!xdr_vector(stream, (char*) data, n,
                   sizeof(float), (xdrproc_t) xdr_float))
   {
      return 1;
   }
   return 0;
}

int xdr_read_double(XDR *stream, double *data, const int n)
{
   if (!xdr_vector(stream, (char*) data, n,
                   sizeof(double), (xdrproc_t) xdr_double))
   {
      return 1;
   }
   return 0;
}

int xdr_write_string(XDR *stream, char **data, int max_length)
{
   if (!xdr_string(stream, data, max_length))
   {
      return 1;
   }
   return 0;
}

// xdr read routines from Scott 
int xdr_write_int(XDR *stream, int *data, const int n)
{
   if (!xdr_vector(stream, (char *) data, n,
                   sizeof(int), (xdrproc_t) xdr_int))
   {
      return 1;
   }
   return 0;
}

int xdr_write_float(XDR *stream, float *data, const int n)
{
   if (!xdr_vector(stream, (char*) data, n,
                   sizeof(float), (xdrproc_t) xdr_float))
   {
      return 1;
   }
   return 0;
}

int xdr_write_double(XDR *stream, double *data, const int n)
{
   if (!xdr_vector(stream, (char*) data, n,
                   sizeof(double), (xdrproc_t) xdr_double))
   {
      return 1;
   }
   return 0;
}



#if 0
void ComputeRange()
{
   int variable;
   int patch;
   
   for(variable = 0; variable < NumVariables; variable ++)
   {
      ScalarRange[variable][0] = VTK_FLOAT_MAX;
      ScalarRange[variable][1] = VTK_FLOAT_MIN;
   }

   for(patch = 0; patch < NumPatches; patch++)
   {
      for(variable = 0; variable < NumVariables; variable ++)
      {
	 float *patchrange;
	 Scalars[patch][variable] -> ComputeRange();
	 
	 patchrange = Scalars[patch][variable] -> GetRange();

	 if ( ScalarRange[variable][0] > patchrange[0]) 
	    ScalarRange[variable][0] = patchrange[0];

	 if ( ScalarRange[variable][1] < patchrange[1]) 
	    ScalarRange[variable][1] = patchrange[1];
      }
   }

   // If no range then create a small artificial one to prevent errors
   // Could be fixed elsewhere but patching this here seemed to be easier
   for(variable = 0; variable < NumVariables; variable ++)
   {
      if ( ScalarRange[variable][0] == ScalarRange[variable][1] )
	    ScalarRange[variable][1] += 0.000001;
   }
}
#endif

int main(int argc, char **argv)
{
   char * FileName;
   int FileFormat;
   
   int NumPatches;
   
   char *VariableName;
   int Variable;
   
   int NumVariables;
   char **VariableNames;
   
   int *PatchLevel;

   int *Size;
   
   int NumLevels;
   
   int Scales[3];
   
   float **ScalarRange;
   
   float Bounds[6];

   float Spacing[3];
   
   float Scaling;

   unsigned int **StartingPosition;

   int NumDims = 3;

   int DataFormat;

   int DataElementSize;
   
   FILE *file;

   FILE *out_file;

   int dim[3];

   int max_level = 0;

   // Should these be of type float?  or double?  both?

   double delta[3];
   double corner[3];

   XDR out_stream;

   int i;
   int ii;

   int patch;

   if ( (file = fopen(argv[1], "r")) == NULL )
   {
      cerr << "Can't open file" << argv[1] << endl;
      perror("Error reported from fopen");
      return 1;
   }

   fscanf(file, "%d %d", &dim[0], &dim[1]);
   dim[2] = 1;

   cout << "Dim = " << dim[0] << "," << dim[1] << endl;

   fscanf(file, "%lf %lf", &delta[0], &delta[1]);

   delta[2] = (delta[0] + delta[1]) / 2.0;

   cout << "Delta = " << delta[0] << "," << delta[1] << endl;

   if ( (out_file = fopen(argv[2], "wb")) == NULL )
   {
      cerr << "Can't open file" << argv[1] << endl;
      perror("Error reported from fopen");
      return 1;
   }

   xdrstdio_create(&out_stream, out_file, XDR_ENCODE);

   // Read in the first int in file; this can be the number of boxes
   // or the file type (if < 0)
   NumPatches=-5;
   xdr_write_int(&out_stream, &NumPatches, 1);

   if ( NumPatches < 0 ) 
   {
      // Dealing with new format 
      FileFormat = abs(NumPatches);
   }
   else
      FileFormat = 0;

   if (FileFormat > 6) 
   {
      cout << "File format is not supported" << endl;
      exit(1);
   }

   switch (FileFormat) 
   {
      case 5:
      {
	 double time;
	 time=0;
	 xdr_write_double(&out_stream, &time, 1);
      }
   }

   
   // SGS

   switch (FileFormat) 
   {
      case 1:
      case 2:
      case 3:
      case 4:
      case 5:
      {
	 // Read in the number of boxes
	 NumPatches=1;
	 xdr_write_int(&out_stream, &NumPatches, 1);
	 break;
      }
   }


   // Default to double for older file formats
   DataFormat=0;
   switch (FileFormat) 
   {
      case 3:
      case 4:
      case 5:
      {
	 // Read storage format for data values
	 DataFormat=1;
	 xdr_write_int(&out_stream, &DataFormat, 1);
	 break;
      }
   }

   switch(DataFormat)
   {
      case 0:
      {
	 DataElementSize = 8;
	 break;
      }
      case 1:
      {
	 DataElementSize = 4;
	 break;
      }
      
   }

   NumVariables=1;
   xdr_write_int(&out_stream, &NumVariables, 1);

   switch (FileFormat) 
   {
      case 1:
      case 2:
      case 3:
      case 4:
      case 5:
      {
         VariableNames = (char **)calloc(NumVariables, sizeof(char *));
         for ( int i = 0; i < NumVariables; i++)
         {
            // Read in the name for the component
	    VariableNames[i]="Population";
            xdr_write_string(&out_stream, &VariableNames[i], MAX_STRING);
         }
         break;
      }
   }

   switch (FileFormat) {
      case 2:
      case 3:
      case 4:
      case 5:
      {
	 // For now ignore this information
	 int *scale;
	 NumLevels=1;
	 xdr_write_int(&out_stream, &NumLevels, 1);

	 scale = (int *)malloc(sizeof(int)*(NumLevels-1)*3);
	 
	 xdr_write_int(&out_stream, scale, (NumLevels-1)*3);
	 
	 free(scale);
	 break;
      }
   }


   PatchLevel = (int *)calloc(NumPatches, sizeof(int));

   for( i = 0; i < NumDims; i++) 
   {
      Spacing[i] = VTK_FLOAT_MAX;
   }

   switch (FileFormat)
   {
      case 4:
      case 5:
      {
	 int NumberOfPatchBoundaries;
	 int BoundaryLevel;
	 
	 NumberOfPatchBoundaries=1;
	 xdr_write_int(&out_stream, &NumberOfPatchBoundaries, 1);
	 for(patch = 0; patch < NumberOfPatchBoundaries; patch ++)
	 {
	    double d_bounds[6];
	    d_bounds[0] = 0;
	    d_bounds[1] = 0;
	    d_bounds[2] = 0;
	    d_bounds[3] = dim[0]*delta[0];
	    d_bounds[4] = dim[1]*delta[1];
	    d_bounds[5] = dim[2]*delta[2];
	    xdr_write_double(&out_stream, d_bounds, 6);

	    BoundaryLevel = 0;
	    xdr_write_int(&out_stream, &BoundaryLevel, 1);
	 }
	 break;
      }
   }


   for(patch = 0; patch < NumPatches; patch++)
   {
      int lower[3];
      int upper[3];
      int size;

      switch (FileFormat) 
      {
         case 1:
	 case 2:
	 case 3:
	 case 4:
	 case 5:
         {
            // Read in the level number; ignore this for now
	    PatchLevel[patch] = 0;
            xdr_write_int(&out_stream, &PatchLevel[patch], 1);
            max_level = PatchLevel[patch] > max_level ? 
               PatchLevel[patch] : max_level;
            break;
         }
      }

      lower[0] = 0;
      lower[1] = 0;
      lower[2] = 0;

      xdr_write_int(&out_stream, lower, 3);

      upper[0] = dim[0] - 1;
      upper[1] = dim[1] - 1;
      upper[2] = 1;

      xdr_write_int(&out_stream, upper, 3);

      switch (FileFormat) 
      {
         case 1:
         case 2:
         case 3:
         case 4:
         case 5:
         {
	    // From input
            xdr_write_double(&out_stream, delta, 3);

	    corner[0] = delta[0] / 2.0;
	    corner[1] = delta[1] / 2.0;
	    corner[2] = delta[2] / 2.0;
	    // SGS mod here
            xdr_write_double(&out_stream, corner, 3);

            break;
         }
      }

      // Number of Cells
      size = dim[0]*dim[1]*dim[2];

      float ar[3];

      // Find the smallest delta spacing in all patches, use that for
      // spacing for the entire dataset
      for( i = 0; i < NumDims; i++) 
      {
         ar[i] = delta[i];

	 if ( ar[i] < Spacing[i] ) 
	    Spacing[i] = ar[i];
      }

      float origin[3];

      switch (DataFormat) 
      {
	 int x, y;
	 float value;
	 // Read doubles
	 case 0:
	 {
	    double *temp = new double[size*NumVariables];
	    xdr_write_double(&out_stream, temp, size*NumVariables);
	    delete[] temp;
	    break;
	 }

	 // Read floats
	 case 1:
	 {
	    float *temp = new float[size*NumVariables];

	    for(x=0; x < dim[0]; x++)
	       for(y=0; y < dim[1]; y++)
	       {
		  fscanf(file, "%f", &value);
		  temp[y*dim[0] + x] = value;
	       }
	    
	    xdr_write_float(&out_stream, temp, size*NumVariables);
	    delete[] temp;
	    break;
	 }
      }
      // End of patch loop
   }

   fflush(out_file);
   xdr_destroy(&out_stream);
   fclose(out_file);

   return 0;
}

