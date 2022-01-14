//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/toolbox/timers/Statistician.h $
// Package:     SAMRAI toolbox
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    \f$       \f$
// Modified:    \f$ \f$
// Description: Singleton manager class for statistic objects.
//

#ifndef included_tbox_Statistician
#define included_tbox_Statistician

#include "SAMRAI_config.h"
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include "tbox/Serializable.h"
#include "tbox/Statistic.h"
#ifndef included_String
#include <string>
#define included_String
#endif

namespace SAMRAI {
   namespace tbox {

class StatisticRestartDatabase;

/**
 * Class Statistician is a Singleton class that manages a simple 
 * database of Statistic objects.  This class provides a single point
 * of access to statistic objects so that any one of them may be updated 
 * or recorded at any point in the code.  Access to the Singleton 
 * statistician instance follows the standard SAMRAI implementation 
 * found in other classes with similar Singleton behavior.  See static
 * member functions below for more information.
 *
 * Statistic objects can be to the database or accessed in code as follows:
 *
 *     Pointer<Statistic> stat = 
 *           Statistician::getStatistician->
 *           getStatistic("name", "PROC_STAT");
 *
 * Here `name' is the name string identifier for the statistic object and 
 * `PROC_STAT' is the type of statistic. See discussion for the getStatistic()
 * method below for more information.
 *
 * The statistic state is saved to restart files when restart file generation 
 * is active.  This allows users to continue to accumulate timing 
 * information when restarting a run.  If desired, all statistics can be 
 * reset when restarting by calling the function resetAllStatistics().
 *
 * A variety of print options exist to dump statistics data.  Notably,
 * the printSpreadSheetOutput("print_dir") will write statistics data in
 * a tab-separated format to files in the supplied directory name.  The
 * naming convention for statistics data is "\<name\>-\<type\>.txt" where \<name\>
 * is the name of the statistic and \<type\> is either proc or patch stat.  The
 * files may be read in to a spreadsheet program such as MS Excel.
 *
 * For more information about data that can be recorded with statistics,
 * consult the header file for the Statistic class.
 *
 * @see tbox::Statistic
 */

class Statistician
{
   friend class StatisticRestartDatabase; 
public:
   /**
    * Create the singleton instance of the statistic manager and return
    * a pointer to it.  This function is provided so that so that 
    * users can control whether statistic information will be written 
    * to/read from restart files.
    *
    * The statistic restart database object is also resistered for writing
    * subsequent restart files when the first boolean argument is true.  
    * Whether the statistic database will write statistics to restart files 
    * during program execution is determined by this argument (true by 
    * default).  Regardless of the value of this argument, statistics that 
    * exist in the restart file will be read from restart when a run is 
    * restarted and the second argument is true.
    *
    * Generally, this routine should only be called once during program
    * execution.  If the statistician has been previously created (e.g.,
    * by an earlier call to this routine) this routine will do nothing
    * other than return the pointer to the existing singleton instance.
    */
   static Statistician* createStatistician(
      bool register_for_restart = true,
      bool read_from_restart = true);

   /**
    * Return a pointer to the singleton statistician instance.
    * All access to the Statistician object is through the
    * getStatistician() function.  For example, to add a statistic 
    * object with the name "my_stat" to the statistician, use the 
    * following call: 
    * Statistician::getStatistician()->addStatistic("my_stat").
    */
   static Statistician* getStatistician();

   /**
    * Deallocate the Statistician instance.  Note that it is not
    * necessary to call freeStatistician() at program termination, since 
    * it is automatically called by the ShutdownRegistry class. 
    */
   static void freeStatistician();

   /**
    * Return pointer to statistic object with the given name string.
    * If a statistics with the given name already exists in the database
    * of statistics, the statistic with that name will be returned.  
    * Otherwise, a new statistic will be created with that name.  The
    * stat_type string identifier is only used when a new statistic
    * object must be created.  The two avaible options are processor
    * statistics and patch statistics which are indicated by the strings
    * "PROC_STAT" and "PATCH_STAT", respectively.
    *
    * When assertion checking is active, an assertion will result if wither
    * string is empty.
    */
   virtual Pointer<Statistic> getStatistic(const std::string& name,
                                                     const std::string& stat_type);

   /**
    * Return true if a statistic whose name matches the argument string
    * exists in the database of statistics controlled by the statistician.
    * If a match is found, the statistic pointer in the argument list is set
    * to that statistic.  Otherwise, return false and return a null pointer.
    * If the name string is empty, a null pointer is returned.
    */
   virtual bool checkStatisticExists(Pointer<Statistic>& stat,
                                     const std::string& name) const;

   /**
    * Return integer number of local processor statistics maintained 
    * by statistician.
    */
   virtual int getNumberProcessorStats() const;

   /**
    * Return integer number of local patch statistics maintained 
    * by statistician.
    */
   virtual int getNumberPatchStats() const;

   /**
    * Reset all processor statistics to contain no information. The primary 
    * intent of this function is to avoid using restarted statistic values 
    * when performing a restarted run.  However, it can be called anytime.
    */
   virtual void resetProcessorStatistics();  

   /**
    * Reset all patch statistics to contain no information. The primary 
    * intent of this function is to avoid using restarted statistic values 
    * when performing a restarted run.  However, it can be called anytime.
    */
   virtual void resetPatchStatistics(); 

   /**
    * Reset all patch and processor statistics to contain no information.
    */
   virtual void resetStatistics(); 

   /**
    * Return integer instance identifier for processor statistic with given
    * name.  If this statistician object maintains no processor statistic
    * with that name, then a warning message results and the return value
    * will be the invalid instance identifier "-1".
    */
   virtual int getProcStatId(const std::string& name) const;

   /**
    * Return number of sequence entries for processor statistic with given 
    * integer identifier.   For convenience, the routine getProcStatId() is 
    * provided to map the statistic string name to the proper integer 
    * identifier.
    *
    * When assertion checking is active, the identifier must be valid or 
    * an assertion will result.
    */
   virtual int getGlobalProcStatSequenceLength(int proc_stat_id);

   /**
    * Return statistic data value for processor statistic with given integer
    * identifier, sequence number, and processor number.   For convenience,
    * the routine getProcStatId() is provided to map the statistic string 
    * name to the proper integer identifier.  The function 
    * getGlobalProcStatSequenceLength() provides the sequence length for
    * a given processor statistic. 
    *
    * When assertion checking is active, the identifier, sequence number, 
    * and processor number must be valid or an assertion will result.
    */
   virtual double getGlobalProcStatValue(int proc_stat_id,
                                         int seq_num, 
                                         int proc_num);

   /**
    * Return global sum of processor statistic with given integer
    * identifier and sequence number.   To identify the correct integer
    * identifier and valid sequence numbers, the method getProcStatId() maps
    * the statistic string name to its integer identifier and the method
    * getGlobalProcStatSequenceLength() returns the maximum sequence length
    * for the processor statistic. 
    *
    * When assertion checking is active, the identifier and sequence number
    * must be valid or an assertion will result.
    */
   virtual double getGlobalProcStatSum(int proc_stat_id,
                                       int seq_num);

   /**
    * Return global max of processor statistic with given integer
    * identifier and sequence number.   To identify the correct integer
    * identifier and valid sequence numbers, the method getProcStatId() maps
    * the statistic string name to its integer identifier and the method
    * getGlobalProcStatSequenceLength() returns the maximum sequence length
    * for the processor statistic. 
    *
    * When assertion checking is active, the identifier and sequence number
    * must be valid or an assertion will result.
    */
   virtual double getGlobalProcStatMax(int proc_stat_id,
                                       int seq_num);

   /**
    * Returns rank of processor holding global max for the processor
    * statistic specified by the given integer identifyer and sequence 
    * number.
    *
    * When assertion checking is active, the identifier and sequence number
    * must be valid or an assertion will result.
    */
   virtual int getGlobalProcStatMaxProcessorId(int proc_stat_id,
                                               int seq_num);

   /**
    * Return global min of processor statistic with given integer
    * identifier and sequence number.   To identify the correct integer
    * identifier and valid sequence numbers, the method getProcStatId() maps
    * the statistic string name to its integer identifier and the method
    * getGlobalProcStatSequenceLength() returns the maximum sequence length
    * for the processor statistic. 
    *
    * When assertion checking is active, the identifier and sequence number
    * must be valid or an assertion will result.
    */
   virtual double getGlobalProcStatMin(int proc_stat_id,
                                       int seq_num);

   /**
    * Returns rank of processor holding global max for the processor
    * statistic specified by the given integer identifyer and sequence 
    * number.
    *
    * When assertion checking is active, the identifier and sequence number
    * must be valid or an assertion will result.
    */
   virtual int getGlobalProcStatMinProcessorId(int proc_stat_id,
                                               int seq_num);

   /**
    * Print global processor statistic data for a particular statistic
    * to given output stream.  Floating point precision may be specified
    * (default is 12).  Note that this method generates a general dump of 
    * the data but does NOT generate it in tabulated form.  To generate 
    * tabulated data, see the printGlobalProcStatDataFormatted() method.
    */
   virtual void printGlobalProcStatData(int proc_stat_id,
                                        std::ostream& os,
                                        int precision = 12); 

   /**
    * Print processor stat data in formatted output to given output 
    * stream.  Floating point precision may be specified (default is 12).  
    */
   virtual void printGlobalProcStatDataFormatted(int proc_stat_id,
                                                 std::ostream& os,
                                                 int precision = 12);

   /**
    * Print stat data for specified processor in formatted output to 
    * given output stream.  Floating point precision may be specified 
    * (default is 12).  
    */
   virtual void printGlobalProcStatDataFormatted(int proc_stat_id,
                                                 int proc_id,
                                                 std::ostream& os,
                                                 int precision = 12);

   /**
    * Return integer instance identifier for patch statistic with given
    * name.  If this statistician object maintains no patch statistic
    * with that name, then a warning message results and the return value
    * will be the invalid instance identifier "-1".
    */
   virtual int getPatchStatId(const std::string& name) const;

   /**
    * Return number of sequence entries for patch statistic with given
    * integer identifier.   For convenience, the routine getPatchStatId() 
    * is provided to map the statistic string name to the proper integer
    * identifier.
    *
    * When assertion checking is active, the identifier must be valid or
    * an assertion will result.
    */
   virtual int getGlobalPatchStatSequenceLength(int patch_stat_id);
  
   /**
    * Return number of patch entries for patch statistic with given
    * integer identifier, and sequence number.   For convenience, the 
    * routine getPatchStatId() is provided to map the statistic string 
    * name to the proper integer identifier.  The function 
    * getGlobalPatchStatSequenceLength() provides the sequence length for
    * a given patch statistic.
    *
    * When assertion checking is active, the identifier and sequence number
    * must be valid or an assertion will result.
    */
   virtual int getGlobalPatchStatNumberPatches(int patch_stat_id,
                                               int seq_num);

   /**
    * Return global processor mapping for patch statistic with given integer
    * identifier, sequence number, and patch number.  For convenience,
    * the routine getPatchStatId() is provided to map the statistic string
    * name to the proper integer identifier.  The function
    * getGlobalPatchStatSequenceLength() provides the sequence length for
    * a given patch statistic.  The function
    * getGlobalPatchStatNumberPatches() gives the number of patches
    * associated with a patch statistic and sequence number.
    *
    * When assertion checking is active, the identifier, sequence number,
    * and patch number must be valid or an assertion will result.
    */
   virtual int getGlobalPatchStatPatchMapping(int patch_stat_id,
                                              int seq_num,
                                              int patch_num);

   /**
    * Return statistic data value for patch statistic with given integer
    * identifier, sequence number, and patch number.   For convenience,
    * the routine getPatchStatId() is provided to map the statistic string
    * name to the proper integer identifier.  The function 
    * getGlobalPatchStatSequenceLength() provides the sequence length for
    * a given patch statistic.  The function 
    * getGlobalPatchStatNumberPatches() gives the number of patches
    * associated with a patch statistic and sequence number.
    *
    * When assertion checking is active, the identifier, sequence number,
    * and patch number must be valid or an assertion will result.
    */
   virtual double getGlobalPatchStatValue(int patch_stat_id,
                                          int seq_num,
                                          int patch_num);

   /**
    * Return global sum of patch statistic with given integer
    * identifier and sequence number.   To identify the correct integer
    * identifier and valid sequence numbers, the method getPatchStatId() maps
    * the statistic string name to its integer identifier and the method
    * getGlobalPatchStatSequenceLength() returns the maximum sequence length
    * for the processor statistic. 
    *
    * When assertion checking is active, the identifier and sequence number
    * must be valid or an assertion will result.
    */
   virtual double getGlobalPatchStatSum(int patch_stat_id,
                                        int seq_num);

   /**
    * Return global max of patch statistic with given integer
    * identifier and sequence number.   To identify the correct integer
    * identifier and valid sequence numbers, the method getPatchStatId() maps
    * the statistic string name to its integer identifier and the method
    * getGlobalPatchStatSequenceLength() returns the maximum sequence length
    * for the processor statistic. 
    *
    * When assertion checking is active, the identifier and sequence number
    * must be valid or an assertion will result.
    */
   virtual double getGlobalPatchStatMax(int patch_stat_id,
                                        int seq_num);

   /**
    * Returns ID of patch holding global max for the patch
    * statistic specified by the given integer identifyer and sequence 
    * number.
    *
    * When assertion checking is active, the identifier and sequence number
    * must be valid or an assertion will result.
    */
   virtual int getGlobalPatchStatMaxPatchId(int patch_stat_id,
                                            int seq_num);

   /**
    * Return global min of patch statistic with given integer
    * identifier and sequence number.   To identify the correct integer
    * identifier and valid sequence numbers, the method getPatchStatId() maps
    * the statistic string name to its integer identifier and the method
    * getGlobalPatchStatSequenceLength() returns the maximum sequence length
    * for the processor statistic. 
    *
    * When assertion checking is active, the identifier and sequence number
    * must be valid or an assertion will result.
    */
   virtual double getGlobalPatchStatMin(int patch_stat_id,
                                        int seq_num);

   /**
    * Returns patch ID of patch holding global min for the patch
    * statistic specified by the given integer identifyer and sequence 
    * number.
    *
    * When assertion checking is active, the identifier and sequence number
    * must be valid or an assertion will result.
    */
   virtual int getGlobalPatchStatMinPatchId(int patch_stat_id,
                                            int seq_num);

   /**
    * Returns the sum of patch statistic information for a particular
    * processor.  The patch statistic is specified by its integer identifyer 
    * and sequence number.
    *
    * When assertion checking is active, the identifier, processor id,
    * and sequence number must be valid or an assertion will result.
    */
   virtual double getGlobalPatchStatProcessorSum(int patch_stat_id,
                                                int processor_id,
                                                int seq_num);

   /**
    * Returns the maximum value of the patch statistic data summed
    * on each processor.  That is, patch statistic information is 
    * summed on each processor, and this method returns the maximum 
    * value, across all processors, of this summed data. The patch 
    * statistic is specified by its integer identifyer and sequence 
    * number.
    *
    * When assertion checking is active, the identifier and sequence 
    * number must be valid or an assertion will result.
    */
   virtual double getGlobalPatchStatProcessorSumMax(int patch_stat_id,
                                                    int seq_num);

   /**
    * Returns the processor ID which holds the maximum value of 
    * summed patch statistic information across processors.  See
    * the discussion for the method getGlobalPatchStatProcessorSumMax()
    * for more information on the summed patch statistic information
    * on processors. 
    *
    * When assertion checking is active, the identifier and sequence 
    * number must be valid or an assertion will result.
    */
   virtual int getGlobalPatchStatProcessorSumMaxId(int patch_stat_id,
                                                   int seq_num);
   /**
    * Returns the minimum value of the patch statistic data summed
    * on each processor.  That is, patch statistic information is 
    * summed on each processor, and this method returns the minimum 
    * value, across all processors, of this summed data. The patch 
    * statistic is specified by its integer identifyer and sequence 
    * number.
    *
    * When assertion checking is active, the identifier and sequence 
    * number must be valid or an assertion will result.
    */
   virtual double getGlobalPatchStatProcessorSumMin(int patch_stat_id,
                                                    int seq_num);

   /**
    * Returns the processor ID which holds the minimum value of 
    * summed patch statistic information across processors.  See
    * the discussion for the method getGlobalPatchStatProcessorSumMin()
    * for more information on the summed patch statistic information
    * on processors. 
    *
    * When assertion checking is active, the identifier and sequence 
    * number must be valid or an assertion will result.
    */
   virtual int getGlobalPatchStatProcessorSumMinId(int patch_stat_id,
                                                   int seq_num);

   /**
    * Return number of patches on the specified processor number for 
    * patch statistic with given identifier, and sequence number.   
    *
    * When assertion checking is active, the identifier and sequence number
    * must be valid or an assertion will result.
    */
   virtual int getGlobalPatchStatNumberPatchesOnProc(int patch_stat_id,
                                                     int seq_num,
                                                     int proc_id);

   /**
    * Returns the maximum number of patches per processor for the
    * specified patch statistic.  
    *
    * When assertion checking is active, the identifier and sequence 
    * number must be valid or an assertion will result.
    */
   virtual int getGlobalPatchStatMaxPatchesPerProc(int patch_stat_id,
                                                   int seq_num);

   /**
    * Returns the processor ID holding the maximum number of patches 
    * per processor for the specified patch statistic.  
    *
    * When assertion checking is active, the identifier and sequence 
    * number must be valid or an assertion will result.
    */
   virtual int getGlobalPatchStatMaxPatchesPerProcId(int patch_stat_id,
                                                     int seq_num);

   /**
    * Returns the minimum number of patches per processor for the
    * specified patch statistic.  
    *
    * When assertion checking is active, the identifier and sequence 
    * number must be valid or an assertion will result.
    */
   virtual int getGlobalPatchStatMinPatchesPerProc(int patch_stat_id,
                                                   int seq_num);

   /**
    * Returns the processor ID holding the minimum number of patches 
    * per processor for the specified patch statistic.  
    *
    * When assertion checking is active, the identifier and sequence 
    * number must be valid or an assertion will result.
    */
   virtual int getGlobalPatchStatMinPatchesPerProcId(int patch_stat_id,
                                                     int seq_num);
   /**
    * Print global processor statistic data for a particular statistic
    * to given output stream.  Floating point precision may be specified
    * (default is 12).  Note that this method generates a general dump of 
    * the data but does NOT generate it in tabulated form.  To generate 
    * tabulated data, see the printGlobalPatchStatDataFormatted() method.
    */
   virtual void printGlobalPatchStatData(int patch_stat_id,
                                         std::ostream& os,
                                         int precision = 12);

   /**
    * Print patch stat data in formatted output to given output 
    * stream.  Floating point precision may be specified (default is 12).
    */
   virtual void printGlobalPatchStatDataFormatted(int patch_stat_id,
                                                  std::ostream& os,
                                                  int precision = 12);

   /**
    * Gather all statistic information to the statistician object 
    * on processor zero.  Typically, this routine is called after 
    * a calculation has completed so that statistic data can be 
    * retrieved, analized, printed to a file, etc.  It is not essential
    * that this routine be called, however, as each "get" and "print"
    * routine checks to see if statistic data has been finalized before
    * it peforms its function.
    */
   virtual void finalize();

   /**
    * Print data to given output stream for local statistics managed 
    * by this statistician object.  Note that no fancy formatting is done.
    * Floating point precision can be specified (default is 12).
    */
   virtual void printLocalStatData(std::ostream& os,
                                   int precision = 12) const;

   /**
    * Print global statistic data information to given output stream.
    * The data will NOT be in tabulated form.  Floating point precision 
    * can be specified (default is 12).
    */
   virtual void printAllGlobalStatData(std::ostream& os,
                                       int precision = 12);

   /**
    * Print sums of all global statistic data information to given 
    * output stream. Floating point precision can be specified (default is 12).
    */
   virtual void printAllSummedGlobalStatData(std::ostream& os,
                                             int precision = 12);

   /**
    * Print sums of all global statistic data information to specified
    * filename. Floating point precision can be specified (default is 12).
    */
   virtual void printAllSummedGlobalStatData(const std::string& filename,
                                             int precision = 12);


   /**
    * Write all statistics data in tab-separated format to files in the 
    * supplied directory name.  The naming convention used is "\<name\>-\<type\>.txt" 
    * where \<name\> is the name of the statistic and \<type\> is either proc or 
    * patch stat.  Floating point precision may be specified (default is 12).
    * The files may be read in to a spreadsheet program such as MS Excel.  If 
    * no directory name is supplied, the files will be written to the directory 
    * where the application is run. 
    */
  virtual void printSpreadSheetOutput(const std::string& dirname = std::string(),
                                      int precision = 12);

   /**
    * Write tab-separated statistics data for specified processor.  This
    * method is identical to "printSpreadSheetOutput()" (above), but only
    * prints data for a single processor.  This may be useful for information
    * that is the same across all processors.  This method will only print
    * processor stats. Any patch stats will be ignored. 
    */
  virtual void printSpreadSheetOutputForProcessor(
     const int proc_id,
     const std::string& dirname = std::string(),
     int precision = 12);

protected:
   /**
    * The constructor for Statistician is protected.  Consistent
    * with the definition of a Singleton class, only a statistician object
    * can have access to the constructor for the class.
    */
   Statistician();

   /**
    * Statistician is a Singleton class; its destructor is protected.
    */
   virtual ~Statistician();

   /**
    * Initialize Singleton instance with instance of subclass.  This function
    * is used to make the singleton object unique when inheriting from this
    * base class.
    */
   void registerSingletonSubclassInstance(
      Statistician* subclass_instance);

   /**
    * During finalize() check statistic information on all processors
    * for consistency before generating arrays of data.
    */
   virtual void checkStatsForConsistency(Array<int>& total_patches);

   /**
    * Return true if a processor statistic whose name matches the 
    * argument string exists in the database of statistics controlled 
    * by the statistician.  If a match is found, the statistic pointer 
    * in the argument list is set to that statistic.  Otherwise, return 
    * false and return a null pointer.
    */
   virtual bool checkProcStatExists(Pointer<Statistic>& stat,
                                    const std::string& name) const;

   /**
    * Return true if a patch statistic whose name matches the
    * argument string exists in the database of statistics controlled
    * by the statistician.  If a match is found, the statistic pointer
    * in the argument list is set to that statistic.  Otherwise, return
    * false and return a null pointer.
    */
   virtual bool checkPatchStatExists(Pointer<Statistic>& stat,
                                     const std::string& name) const;

private:

   /*
    * Gets the current maximum number of statistics.
    *
    * If trying to use more statistics than this value 
    * the arrays should be resized.
    */
   int getMaximumNumberOfStatistics();

   /*
    * Set the maximum number of statistics.
    *
    * This will grow the internal arrays used to store values.
    */
   void setMaximumNumberOfStatistics(const int size);

   /**
    * Static data members to manage the singleton statistician instance.
    */
   static Statistician* s_statistician_instance;
   static bool s_registered_callback;

   static void makeStatisticianInstance(bool register_for_restart = true,
                                        bool read_from_restart = true);

   /*
    * Create and initialize state of restart database.
    */
   void initRestartDatabase(bool register_for_restart,
                            bool read_from_restart);

   /*
    * Internal database class for statistician restart capabilities.  See
    * class declaration below.
    */
   StatisticRestartDatabase* d_restart_database_instance;

   /*
    * Count of statistics registered with the statistician and arrays of 
    * pointers to those statistics.
    */
   int d_num_proc_stats;
   Array< Pointer<Statistic> > d_proc_statistics; 
   int d_num_patch_stats;
   Array< Pointer<Statistic> > d_patch_statistics;

   /*
    * Arrays of global statistic data assembled by the finalize() function.
    *
    * Global processor stat data is assembled as 
    *    array(stat id, seq id, proc id) = proc stat value.
    *
    * Global patch stat data is assembled as 
    *    array(stat_id, seq id, global patch id) = patch stat value.
    * 
    * Global patch stat processor data is assembled as
    *    array(stat_id, seq id, global proc id) = patch stats summed on
    *    different processors.
    *
    * The map of patches to processors is assembled as
    *    array(stat_id, seq id, global patch id) = proc number.
    */
   bool d_must_call_finalize;

   Array< Array< Array<double> > > d_global_proc_stat_data;

   Array< Array< Array<double> > > d_global_patch_stat_data;
   Array< Array< Array<double> > > 
      d_global_patch_stat_proc_data;
   Array< Array< Array<int> > > d_global_patch_stat_mapping;


   /*
    * Internal value used to set and grow arrays for storing
    * statistics.
    */
   static const int DEFAULT_NUMBER_OF_TIMERS_INCREMENT = 128;
};

/*
 * Class StatisticRestartDatabase is a separate class used by the 
 * statistician to provide restart capabilities. Each restartable
 * class must be derived from Serializable.  Since Statistician
 * is a Singleton, its destructor is protected.  StatisticRestartDatabase
 * has a publically accessible destructor.  To avoid improper use of this
 * class, it is privately derived from Serializable, and make it a
 * friend of Statistician.  In this way, its methods can only
 * be accessed via Serializable and Statistician.
 */

class StatisticRestartDatabase : private Serializable
{
   friend class Statistician;
public:
   /*
    * The StatisticRestartDatabase constructor caches a copy of the
    * database object name and registers the object with the restart
    * manager for subsequent restart files if the first boolean argument 
    * is true.  If the run is started froma restart file and the second
    * boolean argument is true, we initialize the statistics from restart.
    */
   StatisticRestartDatabase(const std::string& object_name,
                                 bool register_for_restart,
                                 bool read_from_restart);

   /*
    * The destructor for StatisticRestartDatabase unregisters
    * the database object with the restart manager when so registered.
    */
   virtual ~StatisticRestartDatabase();

   /*
    * Put all statistics and their state in the given restart database.
    * This function is overloaded from Serializable.
    */
   void putToDatabase(Pointer<Database> db);

   /*
    * Construct those statistics saved in the restart database.
    */
   void getFromRestart();

private:
   std::string d_object_name;
   bool d_registered_for_restart;

};

}
}
#endif
