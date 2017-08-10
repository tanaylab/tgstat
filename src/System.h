#ifndef _SYSTEM_H_INCLUDED_
#define _SYSTEM_H_INCLUDED_

// This function returns the actual number of unique bytes (cannot shared with other processes) used by pid.
// This is done by counting Private Dirty memory mappings.
// 
// Note: the result might significantly differ from what is being reported by "top" or "ps" utilities since they count
//       also the memory that can be shared between different processes.
// 
// Note: this function has high run-time cost since it needs to open and to parse system files.
// 
// !!!! COMPATIBILITY: only Linux, on other platforms returns 0. !!!!
size_t get_unique_mem_usage(pid_t pid = 0);

#endif
