//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/toolbox/timers/Timer.I $
// Package:     SAMRAI toolbox
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2286 $
// Modified:    $LastChangedDate: 2008-07-09 09:02:19 -0700 (Wed, 09 Jul 2008) $
// Description: Timing wrapper for various routines.
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace SAMRAI {
   namespace tbox {


/*
 *************************************************************************** 
 *                                                                         *
 * Exclusive start and stop routines for timers.                           *
 *                                                                         *
 * For wallclock time: If we have MPI, we use MPI_Wtime to set the         *
 *                     start/stop point.  If we don't have MPI but do      *
 *                     have access to timer utilities in sys/times.h,      *
 *                     we use the time() utility to set the start/start    *
 *                     point.  If we have neither, we set the wallclock    *
 *                     start/stop time to zero.                            *
 *                                                                         *
 * For user time:      If we have access to timer utilities in sys/times.h,*
 *                     we use the times() utility to compute user and      *
 *                     system start/stop point (passing in the tms struct).*
 *                     If we don't have these utilities, we simply set the *
 *                     user and start/stop times to zero.                  *
 *                                                                         *
 * The timer manager manipulates the exclusive time information;  i.e.     *
 * when the timer's startExclusive and stopExclusive routines are called.  *
 *                                                                         *
 *************************************************************************** 
 */

inline void Timer::startExclusive()
{
   if (d_is_active) {

      Clock::timestamp(d_user_start_exclusive,
                            d_system_start_exclusive,
                            d_wallclock_start_exclusive);

   }
}

inline void Timer::stopExclusive()
{
   if (d_is_active) {
      Clock::timestamp(d_user_stop_exclusive,
                            d_system_stop_exclusive,
                            d_wallclock_stop_exclusive);

      d_wallclock_exclusive +=
         double(d_wallclock_stop_exclusive - d_wallclock_start_exclusive);
      d_user_exclusive += 
         double(d_user_stop_exclusive - d_user_start_exclusive);
      d_system_exclusive +=
         double(d_system_stop_exclusive - d_system_start_exclusive);
   }
}

/*
 *************************************************************************** 
 *                                                                         *
 * Simple utility routines to manipulate timers.                           *
 *                                                                         *
 ***************************************************************************
 */

inline const std::string &Timer::getName() const
{
   return(d_name);
}

inline int Timer::getIdentifier() const
{
   return(d_identifier);
}

inline bool Timer::isActive() const
{
   return(d_is_active);
}

inline bool Timer::isRunning() const
{
   return(d_is_running);
}

inline double Timer::getTotalUserTime() const
{
   return(d_user_total/Clock::getClockCycle());
}

inline double Timer::getTotalSystemTime() const
{
   return(d_system_total/Clock::getClockCycle());
}

inline double Timer::getTotalWallclockTime() const
{
   double ret_val = d_wallclock_total;
#ifndef HAVE_MPI
   ret_val /= Clock::getClockCycle();
#endif
   return(ret_val);
}

inline double Timer::getMaxWallclockTime() const
{
   double ret_val = d_max_wallclock;
#ifndef HAVE_MPI
   ret_val /= Clock::getClockCycle();
#endif
   return(ret_val);
}

inline double Timer::getExclusiveUserTime() const
{
   return(d_user_exclusive/Clock::getClockCycle());
}

inline double Timer::getExclusiveSystemTime() const
{
   return(d_system_exclusive/Clock::getClockCycle());
}

inline double Timer::getExclusiveWallclockTime() const
{
   double ret_val = d_wallclock_exclusive;
#ifndef HAVE_MPI
   ret_val /= Clock::getClockCycle();
#endif
   return(ret_val);
}

inline void Timer::setConcurrentTimer(const int id)
{
   /*
    * Dynamically resize array if needed.
    */
   if(id >= d_concurrent_timers.getSize()) {
      d_concurrent_timers.resizeArray(d_concurrent_timers.getSize() 
				 + DEFAULT_NUMBER_OF_TIMERS_INCREMENT);
   }
   d_concurrent_timers[id] = true;
}

inline bool Timer::isConcurrentTimer(const int id) const 
{
   /* 
    * If array is is valid for the id return value otherswise 
    * it must be false since it was never set.
    */
   if(id < d_concurrent_timers.getSize()) {
      return d_concurrent_timers[id];
   } else {
      return false;
   }
}

inline int Timer::getNumberAccesses() const
{
   return(d_accesses);
}

/*
 ***************************************************************************
 *                                                                         *
 * Private utility routines to manipulate timers in database.              *
 *                                                                         *
 ***************************************************************************
 */

inline void Timer::setInactive()
{
   d_is_active = false;
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif

}
}

