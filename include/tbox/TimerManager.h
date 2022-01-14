//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/toolbox/timers/TimerManager.h $
// Package:     SAMRAI toolbox
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2286 $
// Modified:    $LastChangedDate: 2008-07-09 09:02:19 -0700 (Wed, 09 Jul 2008) $
// Description: Singleton timer manager class.
//

#ifndef included_tbox_TimerManager
#define included_tbox_TimerManager

#include "SAMRAI_config.h"

#ifndef included_iostream
#define included_iostream
#include <iostream>
#endif

#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/List.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/Serializable.h"
#ifndef included_String
#include <string>
#define included_String
#endif
#include "tbox/Timer.h"

namespace SAMRAI {
   namespace tbox {

/*!
 * Class TimerManager is a Singleton class that manages a list of
 * timer objects to do performance analysis in SAMRAI library modules 
 * and application codes built with SAMRAI.     
 *
 * Typically, entries in the input file guide timer invocation and 
 * output generation.  Within the source code, timer objects are retrieved 
 * as follows: 
 *
 *    Pointer<Timer> name_timer = 
 *        TimerManager::getManager->getTimer("name");
 *
 * Here `name' is the name string identifier for the timer.  
 *
 * The input data is summarized as follows.  
 *
 * Required input keys: NONE
 * 
 * Optional input keys and defaults:
 * 
 *
 *    - \b    print_exclusive 
 *       (bool) Specifying whether to track and print exclusive times.  
 *       Exclusive times are convenient for identifying time spent inside 
 *       nested routines.  Exclusive time is measured by turning off the 
 *       parent timer when the timer in the nested routine is started.  
 *       It then turns on the parent timer again once the nested timer is 
 *       stopped. This option should be used with some discretion 
 *       because the extra overhead to manage the nested list of timers 
 *       is between four and seven times more expensive than doing simple 
 *       start/stop operations (as is done with the print_total option). 
 *       For this reason, we leave it off by default. [FALSE] 
 *
 *    - \b    print_total 
 *       (bool) Specifies whether to print total (i.e. non-nested) time.  
 *       This is the least expensive way to time parts of the code.  The 
 *       overhead associated with each start/stop sequence is about 
 *       one-half a millionth of a second. [TRUE]
 *
 *    - \b    print_wall 
 *       (bool) Print wallclock time. [TRUE]
 *
 *    - \b    print_user 
 *       (bool) Print user time (measured by system clock). [FALSE]
 *
 *    - \b    print_sys  
 *       (bool) Print system time (measured by system clock). [FALSE]
 *
 *    - \b    print_processor  
 *       (bool) Print times measured on individual processors. [TRUE]
 *
 *    - \b    print_summed  
 *       (bool) Print time summed across all processors. [FALSE]
 *
 *    - \b    print_max  
 *       (bool) Print maximum time spent on any processor, and the 
 *       processor ID that incurred the time. [FALSE]
 *
 *    - \b    print_concurrent  
 *       (bool) Prints the concurrent timer tree, as determined during 
 *       exclusive timing.  That is, for each timer, it prints the list 
 *       of names of the timers in nested routines it calls. [FALSE] 
 *
 *    - \b    print_percentage  
 *       (bool) Prints the percentage of total time with each timer. [TRUE]
 *
 *    - \b    print_timer_overhead 
 *       (bool) Prints some overhead stats associated with the timers. 
 *       Information like the number of times a start/stop sequence was 
 *       called for the timer, and the predicted overhead time associated 
 *       with the timer.  This is a convenient option to occasionally 
 *       check to make sure the timers themselves are not affecting the 
 *       performance of your calculation. [FALSE]  
 *
 *    - \b    print_threshold  
 *       (double) Timers that use up less than (<EM>print_threshold</EM>)
 *       percent of the overall run time are not printed.  This can be 
 *       a convenient option to limit output if you have many timers 
 *       invoked. [0.25]
 *
 *    - \b    timer_list  
 *       (string array) List of timers to be invoked.  The timers can be 
 *       listed individually in <TT>package::class::method</TT> format or 
 *       the entries may contain wildcards to turn on a set of timers in 
 *       a given package or class:
 *       \verbatim
 *          timer_list = "pkg1::*::*", "pkg2::class2::*", ...
 *       \endverbatim
 *
 * 
 *
 * A sample input file entry might look like:
 *
 * \verbatim
 *
 *    print_exclusive      = TRUE
 *    print_timer_overhead = TRUE
 *    timer_list = "algs::HyperbolicLevelIntegrator::advanceLevel()",
 *                 "mesh::GriddingAlgorithm::*",
 *                 "xfer::*::*"
 *
 * \endverbatim
 *       
 * TimerManager expects timer names to be in a certain format to preserve 
 * the wildcard naming capability (i.e. to turn on entire package or class
 * of timers).  See the PDF document in /SAMRAI/docs/userdocs/timing.pdf 
 * for a discussion of how to add timers that maintain this format as well 
 * as a catalog of available timers currently implemented in the library. 
 *
 * Timing recursive function calls will yeild erroneous results and
 * may lead to memory problems.  We recommend {\em not to use timers 
 * to time recursive function calls}.
 * 
 * @see tbox::Timer
 */

class TimerManager
{
   friend class Timer;
public:

   /*!
    * Create the singleton instance of the timer manager. 
    * If the input database pointer is null, no information will be
    * read from the input file.
    *
    * Generally, this routine should only be called once during program 
    * execution.  If the timer manager has been previously created (e.g.,
    * by an earlier call to this routine) this routine will do nothing.
    */
   static void createManager(Pointer<Database> input_db);

   /*!
    * Return a pointer to the singleton instance of the timer manager.
    * All access to the TimerManager object is through the
    * getManager() function.  For example, to add a timer with the name
    * "my_timer" to the timer manager, use the following call:
    * TimerManager::getManager()->addTimer("my_timer").
    */
   static TimerManager* getManager();

   /*!
    * Deallocate the TimerManager instance. Note that it is not
    * necessary to call freeManager() at program termination, since it is
    * automatically called by the ShutdownRegistry class. 
    */
   static void freeManager();

   /*!
    * Return pointer to timer object with the given name string.
    * If a timer with the given name already appears in the database
    * of timers, the timer with that name will be returned.  Otherwise,
    * a new timer will be created with that name.  Typically, only
    * a timer specified (to be turned on) in the input file will be active.
    * Timer names are parsed according to the input data parsing criteria 
    * (described at the top of this class header).  One may override
    * this functionality by adding the timer with a `true' argument.
    * This argument allows one to override the input file criteria and 
    * turn the timer on anyway.
    *
    * When assertion checking is active, an assertion will result if the 
    * string is empty.
    */
   virtual Pointer<Timer> getTimer(const std::string& name,
                                             bool ignore_timer_input = false);

   /*!
    * Return true if a timer whose name matches the argument string
    * exists in the database of timers controlled by the manager.  If
    * a match is found, the timer pointer in the argument list is set
    * to that timer.  Otherwise, return false and return a null pointer.  
    * If the name string is empty, a null pointer is returned.
    */
   virtual bool checkTimerExists(Pointer<Timer>& timer,
                                 const std::string& name) const;

   /*!
    * Return true if a timer whose name matches the argument string
    * exists in the database of timers and is currently running.
    * Otherwise, return false.
    */
   virtual bool checkTimerRunning(const std::string& name) const;

   /*!
    * Reset the times in all timers to zero.
    */
   virtual void resetAllTimers();

   /*!
    * Print the timing statistics to the specified output stream. 
    */
   virtual void print(std::ostream& os = plog);

protected:
   /*!
    * The constructor for TimerManager is protected.  Consistent
    * with the definition of a Singleton class, only a timer manager object
    * can have access to the constructor for the class.
    */
   TimerManager(Pointer<Database> input_db);

   /*!
    * TimerManager is a Singleton class; its destructor is protected.
    */
   virtual ~TimerManager();

   /*!
    * Initialize Singleton instance with instance of subclass.  This function
    * is used to make the singleton object unique when inheriting from this
    * base class.
    */
   void registerSingletonSubclassInstance(
      TimerManager* subclass_instance);

   /*!
    * Mark given timer as running in timer database.
    * If exclusive time option is set, start exclusive time for given
    * timer. Also stop exclusive time for timer on top of exclusive timer
    * stack and push given timer on to that stack.
    *
    * When assertion checking is active, the timer pointer must be non-null.
    */
   virtual void startTime(Timer* timer);

   /*!
    * Mark given timer as not running in timer database.
    * If exclusive time option is set, stop exclusive time for given timer.
    * Also, pop timer off top of exclusive timer stack and start exclusive
    * timer for new top of stack timer.
    *
    * When assertion checking is active, the timer pointer must be non-null.
    */
   virtual void stopTime(Timer* timer);

private:
   /*
    * Static data members to manage the singleton timer manager instance.
    */
   static TimerManager* s_timer_manager_instance;
   static bool s_registered_callback;

   /*
    * Add timer to either the active or inactive timer array.
    */
   bool checkTimerExistsInArray(Pointer<Timer>& timer,
                                const std::string& name,
                                const Array< Pointer<Timer> >& 
                                timer_array,
                                int array_size) const; 

   /*
    * Print a table of values, using values specified in the timer_values
    * array.  column_titles and column_ids specify which columns in 
    * timer_values to print.
    */
   void printTable(const std::string& table_title,
                   const Array<std::string> column_titles,
                   const Array<std::string> timer_names,
                   const int column_ids[],
                   const double timer_values[][18],
                   std::ostream& os);

   /*
    * Same as above, but also print the int max_processor_id integer array
    * in addition to the double values specified in timer_values.  This 
    * function is used when printing maximum values, in which it is 
    * desirable to print the ID of the processor holding the maximum value.
    */
   void printTable(const std::string& table_title,
                   const Array<std::string> column_titles,
                   const Array<std::string> timer_names,
                   const int max_processor_id[][2],
                   const int max_array_id,
                   const int column_ids[],
                   const double timer_values[][18],
                   std::ostream& os);

   /*
    * Output overhead stats for Timers.  
    */
   void printOverhead(const Array<std::string> timer_names,
                      const double timer_values[][18],
                      std::ostream& os);

   /*
    * Output concurrent tree of Timers.  
    */
   void printConcurrent(std::ostream& os);

   /*
    * Build the timer_names, timer_values, and max_processor_id arrays.
    */
   void buildTimerArrays(
                   double timer_values[][18],
                   int max_processor_id[][2],
                   Array<std::string> timer_names);

   /*
    * Build an ordered list array, organizing timers largest to smallest.
    */
   void buildOrderedList(const double timer_values[][18],
                         const int column,
                         int index[],
                         const int array_size);


   /*
    * Checks timer name to determine if it is specified to be turned
    * on.  If it is, return true.  Otherwise, return false. 
    */
   bool checkTimerInNameLists( const std::string& name );

   /*
    * Evaluate consistency of timer database across processors.
    */
   void checkConsistencyAcrossProcessors();

   /*
    * Private member used by the setupTimerDatabase() function to parse
    * input data for managing timers.    
    */
   void getFromInput(Pointer<Database> input_db);

   /*
    * Private member used by the above routine (processInputStringData)
    * and the addTimer routine to add a timer name to the d_package, 
    * d_class, or d_class_method lists.   
    */
   void addTimerToNameLists( const std::string& name );

   /*
    * Quicksort algorithm specialized for timer array..  This
    * implementation is based off of that provided in "Algorithms in 
    * C++", 3rd Edition, Sedgewick.
    */
   static void quicksort(const Array<double>& a, 
                         int index[], 
                         int lo, int hi);

   /*
    * Simple methods to compute percentages, given two doubles.  
    * Performs check to avoid divide-by-zero cases (i.e. where tot = 0)
    * and caps the percentage at 1000 to avoid output irregularities.
    */
   static int computePercentageInt(const double frac, 
				   const double tot); 

   static double computePercentageDouble(const double frac, 
					 const double tot);



   /*
    * Compute the overhead costs of the timing routines 
    * for active and non-active timers.
    *
    * IMPORTANT:  This is destructive of timers so should only
    * be called in the constructor.
    */
   void computeOverheadConstants(void);
   double computeOverheadConstantActiveOrInactive(bool active);

   /*
    * Clear the registered timers.
    */
   void clearArrays(void);

   /*
    * Gets the current maximum number of timers.
    *
    * If trying to use more times than this value 
    * the arrays should be resized.
    */
   int getMaximumNumberOfTimers();

   /*
    * Set the maximum number of timers.
    *
    * This will grow the internal arrays used to store values.
    */
   void setMaximumNumberOfTimers(const int size);

   /*
    * Static constants used by timer manager.
    */
   static int s_main_timer_identifier;
   static int s_inactive_timer_identifier;

   /*
    * Timer accesss overheads.
    */
   double d_timer_active_access_time;
   double d_timer_inactive_access_time;

   /*
    * Main timer used to time overall run time (time between
    * creation and print, or deletion, of TimerManager.  
    */
   Pointer<Timer> d_main_timer;

   /*
    * Count of timers registered with the timer manager and an 
    * array of pointers to those timers.
    */
   int d_num_timers;
   Array< Pointer<Timer> > d_timers; 

   /*
    * An array of dummy inactive timers is used to record 
    * number of accesses to non-active timers.
    */
   int d_num_inactive_timers;
   Array< Pointer<Timer> > d_inactive_timers;

   /*
    * Timer which measures overall run time.  All other timers are 
    * compared against to report time percentages.  
    */

   /*  
    * Boolean array identifying which timers are turned on.
    */
   tbox::Array<bool> d_running_timers;

   /*
    * The timer manager maintains a list (stack) of running timers so
    * that timer objects can maintain exclusive time information as well
    * total elapsed time.  Total elapsed time is the accumulation of time
    * between start and stop calls for a timer.  Exclusive elapsed time 
    * for a timer is its total elapsed time minus any time spent in other 
    * timers while that timer is running.
    */ 
   List< Timer* > d_exclusive_timer_stack;

   /*
    * Lists of timer names generated from the input database.  These are 
    * used to activate specific timers in the code when a program executes.
    */
   List<std::string> d_package_names;
   List<std::string> d_class_names;
   List<std::string> d_class_method_names;

   /*
    * These values hold the length of the package, class, and class_method
    * lists.  They are stored to avoid multiple calls to getNumberOfItems 
    * to improve efficiency.  
    */
   int d_length_package_names;
   int d_length_class_names;
   int d_length_class_method_names;

   /*
    * Print threshold value.  If a "main" timer is specified, the print
    * routines will cutoff printing any timers with percentage less than
    * this value.  For example, if the threshold value is set to 1.0, any
    * timer reporting less than 1% of wallclock time will not be printed.
    */ 
   double d_print_threshold;

   /*
    * Options used in the print routine:
    * Print exclusive and/or total time. 
    * Defaults:  d_print_exclusive=false, d_print_total=true;
    */
   bool d_print_exclusive;
   bool d_print_total;

   /*
    * Print time on individual processor (to log files), maximum time 
    * across all processors, and summed time across all processors.
    * Defaults:  d_print_processor=true, d_print_max=false, 
    * d_print_summed=false;
    */
   bool d_print_processor;
   bool d_print_max;
   bool d_print_summed;

   /*
    * Print user, system, and wallclock time.
    * Defaults:  d_print_user=false, d_print_sys=false, d_print_wall=true;
    */
   bool d_print_user;
   bool d_print_sys;
   bool d_print_wall;

   /*
    * Print percentages of total, concurrent stats, and overhead stats 
    * (i.e. estimated overhead of timer calls).
    * Defaults:  d_print_percentage=true, d_print_concurrent=false, 
    *            d_print_timer_overhead=true;
    */
   bool d_print_percentage;
   bool d_print_concurrent;
   bool d_print_timer_overhead;

   /*
    * Internal value used to set and grow arrays for storing
    * timers.
    */
   static const int DEFAULT_NUMBER_OF_TIMERS_INCREMENT = 128;
};


}
}
#endif
