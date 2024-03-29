//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/toolbox/inputdb/InputManager.C $
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	An input manager singleton class that parses input files
//

#include "tbox/InputManager.h"
#include <stdlib.h>
#include <stdio.h>
#include "tbox/SAMRAI_MPI.h"
#include "tbox/Parser.h"
#include "tbox/PIO.h"
#include "tbox/SAMRAIManager.h"
#include "tbox/ShutdownRegistry.h"
#include "tbox/Utilities.h"

#ifndef NULL
#define NULL (0)
#endif

namespace SAMRAI {
   namespace tbox {

InputManager *InputManager::s_manager_instance = NULL;
bool InputManager::s_registered_callback = false;

Pointer<Database> InputManager::s_input_db
   = Pointer<Database>(NULL);

/*
*************************************************************************
*									*
* Basic singleton classes to create, set, and destroy the manager	*
* instance.								*
*									*
*************************************************************************
*/

InputManager *InputManager::getManager()
{
   if (!s_manager_instance) {
      s_manager_instance = new InputManager;
   }
   if (!s_registered_callback) {
      ShutdownRegistry::registerShutdownRoutine(freeManager,
				ShutdownRegistry::priorityInputManager);
      s_registered_callback = true;
   }
   return(s_manager_instance);
}

void InputManager::setManager(InputManager *manager)
{
   if (s_manager_instance) {
      delete s_manager_instance;
   }
   if (!s_registered_callback) {
      ShutdownRegistry::registerShutdownRoutine(freeManager,
			     ShutdownRegistry::priorityInputManager);

      s_registered_callback = true;
   }
   s_manager_instance = manager;
}
  
void InputManager::freeManager()
{
   if (s_manager_instance) delete s_manager_instance;
   s_manager_instance = ((InputManager *) NULL);

   s_input_db = ((InputDatabase *) NULL);
}

/*
*************************************************************************
*									*
* The constructor and destructor are protected and call only be called	*
* by the singleton class or its subclasses.				*
*									*
*************************************************************************
*/

InputManager::InputManager()
{
}

InputManager::~InputManager()
{
}

/*
*************************************************************************
*									*
* Return whether or not the manager contains an valid input database.   *
*									*
*************************************************************************
*/

bool InputManager::inputDatabaseExists()
{
   return(!(s_input_db.isNull()));
}

/*
*************************************************************************
*									*
* Parse the specified input file and return the new database.		*
*									*
*************************************************************************
*/

Pointer<InputDatabase>
InputManager::parseInputFile(const std::string& filename)
{
   Pointer<InputDatabase> db = new InputDatabase("main");
   this->parseInputFile(filename, db);
   return(db);
}


/*
*************************************************************************
*									*
* Accessor method for InputManger's root input database.                *
*									*
*************************************************************************
*/
Pointer<Database> InputManager::getInputDatabase() 
{
   return(s_input_db);
}

/*
*************************************************************************
*									*
* Parse the specified input file into the given database.		*
*									*
*************************************************************************
*/

void InputManager::parseInputFile(
   const std::string& filename, Pointer<InputDatabase> db)
{
   FILE* fstream = NULL;
   if (SAMRAI_MPI::getRank() == 0) {
      fstream = fopen(filename.c_str(), "r");
   }
   int worked = (fstream ? 1 : 0);
#ifdef HAVE_MPI
   worked = SAMRAI_MPI::bcast(worked, 0);
#endif
   if (!worked) {
      TBOX_ERROR("InputManager:: Could not open input file``" <<
                  filename.c_str() << "''\n");
   }

   /*
    * Parse input file.
    */
   Parser *parser = new Parser();
   const int errors = parser->parse(filename, fstream, db);
   const int warnings = parser->getNumberWarnings();

   if (errors > 0) {
      TBOX_WARNING("InputManager: Errors = " << errors
                   << ", Warnings = " << warnings
                   << "\n when parsing input file = " << filename << std::endl);
      db->printClassData(plog);
      TBOX_ERROR("InputManager exiting..." << std::endl);
   }
   if (warnings > 0) {
      TBOX_WARNING("InputManager: Warnings  = " << warnings
                   << "\n when parsing input file = " << filename << std::endl);
   }

   /*
    * Store the root database in the static s_input_db variable.
    */
   s_input_db = db;

   delete parser;
   if (fstream) fclose(fstream);
}

}
}
