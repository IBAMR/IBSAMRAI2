//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/tools/restart/main.C $
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1965 $
// Modified:    $LastChangedDate: 2008-02-07 12:41:09 -0800 (Thu, 07 Feb 2008) $
// Description: Main program restart-redistribute tool.
//

#include "SAMRAI_config.h"



#include "RedistributedRestartUtility.h"
#include "tbox/HDFDatabase.h"
#include "tbox/RestartManager.h"

#ifndef _MSC_VER
#include <dirent.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <unistd.h>
using namespace std;



// Headers for basic SAMRAI objects

#include "tbox/HDFDatabase.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/RestartManager.h"
#include "tbox/SAMRAIManager.h"

#ifdef HAVE_HDF5

#ifndef NAME_BUFSIZE
#define NAME_BUFSIZE (32)
#endif

using namespace SAMRAI;

#include <iostream>
#include <cstring>
using namespace std;

#ifdef _MSC_VER
//#include <config.h>
#include <string.h>
#include <windows.h>
#include <stdlib.h>

struct dirent { char d_name[1]; };

int scandir(const char *dirname, struct dirent ***namelist,
    int (*select)(struct dirent *),
    int (*compar)(struct dirent **, struct dirent **)) {
  int len;
  char *findIn, *d;
  WIN32_FIND_DATA find;
  HANDLE h;
  int nDir = 0, NDir = 0;
  struct dirent **dir = 0, *selectDir;
  unsigned long ret;

  len    = strlen(dirname);
  findIn = (char*)malloc(len+5);
  strcpy(findIn, dirname);
  for (d = findIn; *d; d++) if (*d=='/') *d='\\';
  if ((len==0)) { strcpy(findIn, ".\\*"); }
  if ((len==1)&& (d[-1]=='.')) { strcpy(findIn, ".\\*"); }
  if ((len>0) && (d[-1]=='\\')) { *d++ = '*'; *d = 0; }
  if ((len>1) && (d[-1]=='.') && (d[-2]=='\\')) { d[-1] = '*'; }

  if ((h=FindFirstFile(findIn, &find))==INVALID_HANDLE_VALUE) {
    ret = GetLastError();
    if (ret != ERROR_NO_MORE_FILES) {
      // TODO: return some error code
    }
    *namelist = dir;
    return nDir;
  }
  do {
    selectDir=(struct dirent*)malloc(sizeof(struct dirent)+strlen(find.cFileName));
    strcpy(selectDir->d_name, find.cFileName);
    if (!select || (*select)(selectDir)) {
      if (nDir==NDir) {
	struct dirent **tempDir = (struct dirent **)calloc(sizeof(struct dirent*), NDir+33);
	if (NDir) memcpy(tempDir, dir, sizeof(struct dirent*)*NDir);
	if (dir) free(dir);
	dir = tempDir;
	NDir += 32;
      }
      dir[nDir] = selectDir;
      nDir++;
      dir[nDir] = 0;
    } else {
      free(selectDir);
    }
  } while (FindNextFile(h, &find));
  ret = GetLastError();
  if (ret != ERROR_NO_MORE_FILES) {
    // TODO: return some error code
  }
  FindClose(h);

  free (findIn);

  if (compar) qsort (dir, nDir, sizeof(*dir),
		     (int(*)(const void*, const void*))compar);

  *namelist = dir;
  return nDir;
}

int alphasort (struct dirent **a, struct dirent **b) {
  return strcmp ((*a)->d_name, (*b)->d_name);
}

#endif

int main( int argc, char *argv[])
{

   const string slash = "/";

   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::startup();

   string read_dirname;
   string write_dirname;
   int restore_num = 0;
   int num_output_files = 1;

   if ( (argc != 5) ) {
	 tbox::pout << "USAGE:  " << argv[0] << " input-dir "
	      << "output-dir restore-number num-output-files\n"
	      << endl;
	 exit(-1);
	 return (-1);
   } else {
      read_dirname = argv[1];
      write_dirname = argv[2];
      restore_num = atoi(argv[3]);
      num_output_files = atoi(argv[4]);
   }

   // make string "input-dir/restore.*****
   char restore_buf[NAME_BUFSIZE];

   sprintf(restore_buf,"/restore.%06d",restore_num);

   string restore_dirname;
   if( read_dirname.compare(0, 1, "/") == 0 ) {
      restore_dirname = read_dirname + restore_buf;
   } else {
      restore_dirname = "."  + slash + read_dirname + restore_buf;
   }

   if( write_dirname.compare(0, 1, "/") != 0 ) {
      write_dirname = "." + slash + write_dirname;
   }

   struct dirent **namelist;

   //directory should have three entries:  ., .., and nodes.*****
   int num_entries = scandir(restore_dirname.c_str(), &namelist, 0, 0); 

   if(num_entries < 0) {
      TBOX_ERROR("restore directory not found.");
   }

   // Expect only a single run to be in restore directory
   if (num_entries > 3) {
      TBOX_ERROR("restore directory should contain restart files for a single run; this probably indicates runs with different number of nodes have been done in in the restart directory");
   }

   string nodes_dirname;
   string prefix = "nodes.";
   for(int i = 0; i < num_entries; i++) {
      if ( strncmp(prefix.c_str(), namelist[i] -> d_name, prefix.length()) == 0 ) {
	 nodes_dirname = namelist[i]->d_name; 
      }
   }

   int num_input_files = 0;

   // Check if nodes_dirname is valid and extract number of processors for the saved run.
   if ( nodes_dirname.size() == 11 ) {

      string int_str = &(nodes_dirname.c_str()[6]);

      num_input_files = atoi(int_str.c_str());

   } else {
      TBOX_ERROR("nodes.***** subdirectory not found in restore directory.  A directory with a name such as nodes.00016 must be present, with the number indicating the number of process for the run");
   }


   free(namelist);

   string full_nodes_dirname = restore_dirname + slash + nodes_dirname;
   num_entries = scandir(full_nodes_dirname.c_str(), &namelist, 0, 0);
   if (num_entries != num_input_files+2) {
      TBOX_ERROR("number of files in nodes subdirectory does not match the number indicated in the directory's name"); 
   }

   free(namelist);

   // file_mapping will have size equal to the lesser value of
   // num_input_files and num_output_files.
   tbox::Array< tbox::Array<int> > file_mapping;
   int file_ratio;
   int remainder;
   if (num_output_files > num_input_files) {
      file_mapping.resizeArray(num_input_files);
      file_ratio = num_output_files / num_input_files;
      remainder = num_output_files % num_input_files;
   } else {
      file_mapping.resizeArray(num_output_files);
      file_ratio = num_input_files / num_output_files;
      remainder = num_input_files % num_output_files;
   }

   int file_counter = 0;

   // fill file_mapping.
   int i;
   for (i = 0; i < file_mapping.size(); i++) {
      if (i < remainder) {
         file_mapping[i].resizeArray(file_ratio+1);
         for (int j = 0; j <= file_ratio; j++) {
            file_mapping[i][j] = file_counter + j;
         }
         file_counter += file_ratio + 1;
      } else {
         file_mapping[i].resizeArray(file_ratio);
         for (int j = 0; j < file_ratio; j++) {
            file_mapping[i][j] = file_counter + j;
         }
         file_counter += file_ratio;
      }
   }

   //Write all of the redistributed restart files.
   RedistributedRestartUtility::writeRedistributedRestartFiles(
      write_dirname,
      read_dirname,
      num_input_files,
      num_output_files,
      file_mapping,
      restore_num);

   tbox::SAMRAIManager::shutdown();

   return(0);

}

#else

int main( int argc, char *argv[])
{
   std::cerr << "This utility requires HDF to work, it was not found when compiling" << std::endl;

   return(-1);
}

#endif
