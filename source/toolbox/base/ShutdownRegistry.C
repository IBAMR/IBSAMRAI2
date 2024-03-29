//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/toolbox/base/ShutdownRegistry.C $
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	Registry of shutdown routines to be called at program exit
//

#include "tbox/ShutdownRegistry.h"

namespace SAMRAI {
   namespace tbox {


#ifndef NULL
#define NULL (0)
#endif

struct ShutdownRegistryItem {
   void (*callback)();
   ShutdownRegistryItem* next;
};

bool ShutdownRegistry::s_initialized = false;

ShutdownRegistryItem* 
ShutdownRegistry::s_shutdown_list[s_number_of_priorities];
ShutdownRegistryItem* 
          ShutdownRegistry::s_shutdown_list_last[s_number_of_priorities];
int ShutdownRegistry::s_num_shutdown_items[s_number_of_priorities];

void ShutdownRegistry::registerShutdownRoutine(void (*callback)(), unsigned char priority)
{
   if(!s_initialized) {
      initialize();
   }

   ShutdownRegistryItem* item = new ShutdownRegistryItem;
   item->callback = callback;
   item->next = (ShutdownRegistryItem*)NULL; 
   if (s_num_shutdown_items[priority] == 0) {
      s_shutdown_list[priority] = item;
   } else {
      s_shutdown_list_last[priority]->next = item;
   }
   s_shutdown_list_last[priority] = item;
   s_num_shutdown_items[priority]++;
}

void ShutdownRegistry::callRegisteredShutdowns()
{
   // only shutdown if something was registered
   if(s_initialized) {

      /* 
       * Shutdown routines may end up registering other shutdown 
       * routines so repeatedly search the list until
       * no more shutdown routines are found.
       */
      bool found = true;
      while(found) {
	 found = false;
	 for(int priority = s_number_of_priorities-1; priority > -1; priority--) {
	    while (s_num_shutdown_items[priority] > 0) {
	       ShutdownRegistryItem* item = s_shutdown_list[priority];
	       (*item->callback)();
	       ShutdownRegistryItem* byebye = item;
	       s_shutdown_list[priority] = item->next;
	       s_num_shutdown_items[priority]--;
	       delete byebye;
	       found = true;
	    }
	    s_shutdown_list_last[priority] = (ShutdownRegistryItem*)NULL;
	 }
      }
   }
}

void ShutdownRegistry::initialize()
{
   for(int priority = s_number_of_priorities-1; priority > -1; priority--) {
      s_shutdown_list[priority] = (ShutdownRegistryItem*)NULL;
      s_shutdown_list_last[priority] = (ShutdownRegistryItem*)NULL;
      s_num_shutdown_items[priority] = 0;
   }

   s_initialized = true;
}


}
}
