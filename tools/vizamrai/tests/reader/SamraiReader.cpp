//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/tools/vizamrai/tests/reader/SamraiReader.cpp $
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
   int *VariableDepths;
   
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

   int dim[3];

   int max_level = 0;

   // Should these be of type float?  or double?  both?

   double delta[3];
   double corner[3];

   XDR stream;

   int i;
   int ii;

   int patch;

   // Sample Reading 
   if ( (file = fopen(argv[1], "rb")) == NULL )
   {
      cerr << "Can't open file" << argv[1] << endl;
      perror("Error reported from fopen");
      return 1;
   }

   xdrstdio_create(&stream, file, XDR_DECODE);

   // Read in the first int in file; this can be the number of boxes
   // or the file type (if < 0)
   xdr_read_int(&stream, &NumPatches, 1);

   if ( NumPatches < 0 ) 
   {
      // Dealing with new format 
      FileFormat = abs(NumPatches);
   }
   else
      FileFormat = 0;

   cout << "FileFormat: " << FileFormat << endl;

   if (FileFormat > 6) 
   {
      cout << "File format is not supported" << endl;
      exit(1);
   }

   switch (FileFormat) 
   {
      case 5:
      case 6:
      {
	 double time;
	 xdr_read_double(&stream, &time, 1);
	 cout << "Plot Time: " << time << endl;
      }
   }


   switch (FileFormat) 
   {
      case 1:
      case 2:
      case 3:
      case 4:
      case 5:
      case 6:
      {
	 // Read in the number of boxes
	 xdr_read_int(&stream, &NumPatches, 1);
	 break;
      }
   }


   cout << "NumPatches: " << NumPatches << endl;

   // Default to double for older file formats
   DataFormat=0;
   switch (FileFormat) 
   {
      case 3:
      case 4:
      case 5:
      case 6:
      {
	 // Read storage format for data values
	 xdr_read_int(&stream, &DataFormat, 1);
	 break;
      }
   }

   switch(DataFormat)
   {
      case 0:
      {
	 DataElementSize = 8;
	 cout << "DataFormat: Doubles " << endl;
	 break;
      }
      case 1:
      {
	 DataElementSize = 4;
	 cout << "DataFormat: floats  " << endl;
	 break;
      }
      
   }

   xdr_read_int(&stream, &NumVariables, 1);
   cout << "NumVariables: " << NumVariables << endl;

   switch (FileFormat) 
   {
      case 1:
      case 2:
      case 3:
      case 4:
      case 5:
      case 6:
      {
         VariableNames = (char **)calloc(NumVariables, sizeof(char *));
         VariableDepths = (int *)calloc(NumVariables, sizeof(int));
         for ( int i = 0; i < NumVariables; i++)
         {
	    // Read in the name for the component
	    xdr_read_string(&stream, &VariableNames[i], MAX_STRING);
	    cout << "VariableNames[" << i << "]: " << VariableNames[i] << endl;
	    switch (FileFormat) 
	    {
	       case 1:
	       case 2:
	       case 3:
	       case 4:
	       case 5:
	       {
		  VariableDepths[i] = 1;
		  cout << "VariableDepths[" << i << "]: " << VariableDepths[i] << endl;
		  
		  break;
	       }

	       case 6:
	       {
		  xdr_read_int(&stream, &VariableDepths[i], 1);
		  cout << "VariableDepths[" << i << "]: " << VariableDepths[i] << endl;
		  break;
	       }
	    }
         }
         break;
      }
   }

   switch (FileFormat) {
      case 2:
      case 3:
      case 4:
      case 5:
      case 6:
      {
	 // For now ignore this information
	 int *scale;
	 xdr_read_int(&stream, &NumLevels, 1);

	 cout << "NumLevels: " << NumLevels << endl;
	 
	 scale = (int *)malloc(sizeof(int)*(NumLevels-1)*3);
	 
	 xdr_read_int(&stream, scale, (NumLevels-1)*3);

	 for(ii = 0; ii < NumLevels-1; ii++)
	 {
	    cout << "scale[" << ii << "]: <"  
		 << scale[ii*3] << ","
		 << scale[ii*3+1] << ","
		 << scale[ii*3+2] << ">" << endl;
	 }
	 
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
      case 6:
      {
	 int NumberOfPatchBoundaries;
	 int BoundaryLevel;

	 xdr_read_int(&stream, &NumberOfPatchBoundaries, 1);
	 for(patch = 0; patch < NumberOfPatchBoundaries; patch ++)
	 {
	    double d_bounds[6];
	    xdr_read_double(&stream, d_bounds, 6);
	    cout << "Patch Boundary: "  << patch << endl;
	    for(ii=0; ii < 3; ii++) 
	    {
	       cout << "    Lower[" << ii << "]: " << d_bounds[ii] << endl;
	    }

	    for(ii=3; ii < 6; ii++) 
	    {
	       cout << "    Upper[" << ii << "]: " << d_bounds[ii] << endl;
	    }

	    xdr_read_int(&stream, &BoundaryLevel, 1);
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
	 case 6:
         {
            // Read in the level number; ignore this for now
            xdr_read_int(&stream, &PatchLevel[patch], 1);
	    cout << "PatchLevel[ " << patch << "]: " << PatchLevel[patch] 
		 << endl;
            max_level = PatchLevel[patch] > max_level ? 
               PatchLevel[patch] : max_level;
            break;
         }
      }

      xdr_read_int(&stream, lower, 3);
      for(ii=0; ii < 3; ii++) 
      {
	 cout << "Lower[" << ii << "]: " << lower[ii] << endl;
      }

      xdr_read_int(&stream, upper, 3);
      for(ii=0; ii < 3; ii++) 
      {
	 cout << "Upper[" << ii << "]: " << upper[ii] << endl;
      }

      switch (FileFormat) 
      {
         case 1:
         case 2:
         case 3:
         case 4:
         case 5:
	 case 6:
         {
            xdr_read_double(&stream, delta, 3);
	    for(ii=0; ii < 3; ii++) 
	    {
	       cout << "Delta[" << ii << "]: " << delta[ii] << endl;
	    }

            xdr_read_double(&stream, corner, 3);
	    for(ii=0; ii < 3; ii++) 
	    {
	       cout << "Corner[" << ii << "]:" << corner[ii] << endl;
	    }

            break;
         }
      }
      dim[0] = upper[0]-lower[0]+1;
      dim[1] = upper[1]-lower[1]+1;
      dim[2] = upper[2]-lower[2]+1;


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

      for( i = 0; i < NumDims; i++)
	  {
         origin[i] = lower[i] * delta[i] + corner[i] - delta[i] * 0.5;
		 cout << "Computed Origin[" << i << "]:" << origin[i] << endl;
	  }


      // Skip over data to next patch
      int skip = 0;
      for(i = 0; i < NumVariables; i++)
      {
	 skip = skip + size * VariableDepths[i];
      }
      skip = skip * DataElementSize;
      int position = xdr_getpos(&stream) + skip;
      cout << "Skipping " << skip << " bytes now at " << position
	   << endl;

      xdr_setpos(&stream, position);

      // End of patch loop
   }

   xdr_destroy(&stream);
   fclose(file);

   return 0;
}

