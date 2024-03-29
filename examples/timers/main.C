//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/examples/timers/main.C $
// Package:     SAMRAI application
// Copyright:   (c) 1997-2002 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Example program to demonstrate timers.
//

#include "SAMRAI_config.h"

#include <string>
using namespace std;

// Headers for basic SAMRAI objects used in this code.
#include "tbox/SAMRAIManager.h"
#include "tbox/Database.h"
#include "tbox/InputManager.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"

using namespace SAMRAI;

int main( int argc, char *argv[] )
{
   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::startup();

   {
      tbox::PIO::logAllNodes("Timer.log");      

#ifdef HAVE_TAU
      TAU_PROFILE("main()", "int (int, char **)", TAU_DEFAULT);
      TAU_PROFILE_SET_NODE(tbox::SAMRAI_MPI::getRank());
#endif


      string input_filename;

      if ( (argc != 2) ) {
	 tbox::pout << "USAGE:  " << argv[0] << " <input filename> "
	      << "  options:\n"
	      << "  none at this time"
	      << endl;
	 tbox::SAMRAI_MPI::abort();
	 return (-1);
      } 

      input_filename = argv[1];

      /*
       * Create an input database "input_db" and parse input file (specified
       * on the command line.
       */
      tbox::Pointer<tbox::Database> input_db = new tbox::InputDatabase("input_db");
      tbox::InputManager::getManager()->parseInputFile(input_filename,input_db);
#if (TESTING == 1)
      tbox::SAMRAIManager::getFromInput(input_db);
#endif

      /*
       * Create timer manager, reading timer list from input.
       */
      tbox::TimerManager::createManager(input_db->getDatabase("TimerManager"));

      /*
       * Add a timer.  All timers are maintained by the tbox::TimerManager
       * and each timer is accessed by its name.  A static pointer is used
       * to avoid the name lookup cost each time it is called (its only
       * looked up the first time).  
       */
      string name = "main::test";
      tbox::Pointer<tbox::Timer> timer = tbox::TimerManager::getManager()->
	 getTimer(name);

      /*
       * Start timer.  If the timer name was not specified in the input
       * file, the timer start() and stop() functions will simply drop through
       * without recording the time.
       */
      timer->start();

      /*
       * Sleep for 10 sec.
       */
      sleep(10);

      /*
       * Stop the timer
       */
      timer->stop();

      /*
       * Print results to log file.
       */
      tbox::TimerManager::getManager()->print(tbox::plog);

   }


   /*
    * We're done.  Shut down application ...
    */
   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAI_MPI::finalize();
   return(0);
}


